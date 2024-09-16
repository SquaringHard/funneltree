#define _USE_MATH_DEFINES   // M_PI
#include "ft.h"
#include <stdexcept>    // runtime_error
#include <algorithm>    // find, swap
#include <string>       // to_string
#include <cmath>        // acos, sqrt, INFINITY
#include <fstream>      // ifstream

TriangleMesh::TriangleMesh(const vector<Point> &ps, const vector<Triangle> &ts) : triangles(ts), points(ps), dictVertices(ps.size()) {
    if (ps.size() > MAX_INDEX) throw out_of_range("too many points");
    if (ts.size() > MAX_INDEX) throw out_of_range("too many faces");

    for (vector<Point>::const_iterator i = points.begin(); i != points.end(); i++) if (find(points.begin(), i, *i) != i)
        throw invalid_argument("point " + to_string(i - points.begin()) + " has duplicates");

    const indexType v = ps.size(), f = ts.size();
    dictVertices.reserve(v);
    dictEdges.reserve(v + f - 2);       // Euler's characteristic
    for (indexType i = 0; i < f; i++) { // add edges and triangle indexes to dictEdges
        for (const Edge &e : {Edge(ts[i].a, ts[i].b), Edge(ts[i].b, ts[i].c), Edge(ts[i].c, ts[i].a)}) {
            const pair<DictEdgeType::iterator, bool> temp = dictEdges.try_emplace(e, array<indexType, 2>{i, MAX_INDEX});
            if (temp.second) continue;

            array<indexType, 2> &eFaces = temp.first->second;
            if (eFaces[1] != MAX_INDEX) throw logic_error("faces " + to_string(eFaces[0]) + ' ' + to_string(eFaces[1]) + ' ' + to_string(i) +
                                                          "occupying same edge");

            eFaces[1] = i;
        }

        for (const indexType j : {ts[i].a, ts[i].b, ts[i].c}) dictVertices[j].push_back(i);
    }

    for (const pair<Edge, array<indexType, 2>> &temp : dictEdges) if (temp.second[1] == MAX_INDEX)
        { const Edge e = temp.first; throw logic_error("floating edge " + to_string(e.a) + '-' + to_string(e.b)); }

    for (indexType i = 0; i < v; i++) if (dictVertices[i].empty()) throw logic_error("floating point " + to_string(i));
}

inline double angle(const double ab, const double bc, const double ca) { return acos((ab * ab + bc * bc - ca * ca) / (ab * bc * 2)); }
inline double calPV(const double pqv, const double pq, const double qv) { return sqrt(pq * pq + qv * qv - pq * qv * cos(pqv) * 2); }

double TriangleMesh::pistance(const indexType a, const indexType b) const {
    const double abx = points[a].x - points[b].x, aby = points[a].y - points[b].y, abz = points[a].z - points[b].z;
    return sqrt(abx * abx + aby * aby + abz * abz);
}

double TriangleMesh::pangle(const indexType a, const indexType b, const indexType c) const { return angle(pistance(a, b), pistance(b, c), pistance(c, a)); }

void Funnel::remove() {
    if (!children) { removed = true; return; }
    children->remove();
    (children + 1)->remove();
}

#pragma omp declare reduction(merge: vector<Funnel*>: omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
typedef unordered_map<Triangle, const Funnel*, HashNComp, HashNComp> FunnelDict;
vector<Funnel*> FunnelTree(const TriangleMesh& mesh, const indexType s) {
    const vector<indexType> &facesAt_s = mesh.dictVertices[s];
    const indexType n = facesAt_s.size();
    vector<Funnel*> list;
    list.reserve(n);
    for (indexType i = 0; i < n; i++) {  // if parallelizing this, read deleteFunnelTree() first
        const Triangle pqv = mesh.triangles[facesAt_s[i]];
        indexType p, q;
        if (pqv.a == s) {
            q = pqv.b;
            p = pqv.c;
        } else if (pqv.b == s) {
            q = pqv.c;
            p = pqv.a;
        } else {
            q = pqv.a;
            p = pqv.b;
        }

        const double spq = mesh.pangle(s, p, q), psw = mesh.pangle(p, s, q);

        // insert all faces containing s instead of just 1 face to sequence so that the funnels never reach these faces again
        // because then the funnels don't have to reach the vertices on these faces
        vector<indexType> newSequence = facesAt_s;
        swap(newSequence[i], newSequence.back());
        list.push_back(new Funnel(p, q, p, newSequence, mesh.pistance(s, p), mesh.pistance(p, q), spq, psw, psw));
    }

    FunnelDict twoChildrenFunnels;
    for (size_t curr = 0;;) {
        const size_t end = list.size();
        vector<Funnel*> next_lvl;
        #pragma omp parallel for reduction(merge: next_lvl) schedule(dynamic)
        for (size_t i = curr; i < end; i++) {
            Funnel *const funnel = list[i];
            if (funnel->removed) continue;

            indexType p = funnel->p, q = funnel->q, x = funnel->x, v;
            double sp = funnel->sp, topright_angle = funnel->topright_angle, pq = funnel->pq, spq = funnel->spq, psq = funnel->psq, psw = funnel->psw;
            vector<indexType> &sequence = funnel->sequence;

            find_v_suchthat_funnelhas2children:

            const array<indexType, 2> eFaces = mesh.dictEdges.at(Edge(x, q));
            const indexType nextFace = eFaces[0] == sequence.back() ? eFaces[1] : eFaces[0];
            if (find(sequence.begin(), sequence.end(), nextFace) != sequence.end()) continue;   // move to next i

            sequence.push_back(nextFace);
            const Triangle temp2 = mesh.triangles[nextFace];
            for (const indexType vNew : {temp2.a, temp2.b, temp2.c}) if (vNew != x && vNew != q) { v = vNew; break; }
            
            topright_angle += mesh.pangle(x, q, v);
            short int sign = 1;
            if (topright_angle >= M_PI) { topright_angle = M_PI * 2 - topright_angle; sign = -1; }
            const double vq = mesh.pistance(v, q), pv = calPV(topright_angle, pq, vq), spv = spq + angle(pv, pq, vq) * sign;

            if (spv >= M_PI) { x = v; goto find_v_suchthat_funnelhas2children; }

            const double sv = calPV(spv, sp, pv), psv = angle(sp, sv, pv), pvq = angle(pv, vq, pq), psw_new = min(psw, psv);
            topright_angle = max(mesh.pangle(x, v, q) - pvq * sign, 0.0);

            if (!(psv < psw)) { // bet ur temping to write psv >= psw
                make_only_Fpv:
                    q = v; pq = pv; spq = spv; psq = psv; psw = psw_new;
                    goto find_v_suchthat_funnelhas2children;
            }

            const double pvs = angle(pv, sv, sp), vsq = psq - psv, vsw = psw - psv, svq = pvq - pvs;
            FunnelDict::const_iterator temp;
            #pragma omp critical (access_pvq)
            temp = twoChildrenFunnels.find({p, v, q});
            if (temp != twoChildrenFunnels.end()) { // clip_off_funnel
                const Funnel *const oldFunnel = temp->second;
                Funnel *const oldChild0 = oldFunnel->children;
                const double sv2 = (oldChild0 + 1)->sp, pvs2 = oldFunnel->pvs;

                if (sv2 > sv) {
                    if (pvs > pvs2) (oldChild0 + 1)->remove();
                    else oldChild0->remove();
                } else {
                    if (sv > sv2) { if (!(pvs > pvs2)) goto make_only_Fpv; }
                    else if (pvs > pvs2) (oldChild0 + 1)->remove();
                    else { oldChild0->remove(); goto make_only_Fpv; }

                    // make only Fvq
                    p = v; x = v; sp = sv; pq = vq; spq = svq; psq = vsq; psw = vsw; topright_angle = 0;
                    goto find_v_suchthat_funnelhas2children;
                }
            }

            // make Fpv and Fvq:
            funnel->pvs = pvs;
            funnel->children = new Funnel[2]{Funnel(p, v, x, sequence, sp, pv, spv, psv, psw_new, topright_angle),
                                             Funnel(v, q, v, sequence, sv, vq, svq, vsq, vsw)};
            #pragma omp critical (access_pvq)
            twoChildrenFunnels[{p, v, q}] = funnel;
            next_lvl.push_back(funnel->children);
            next_lvl.push_back(funnel->children + 1);
        }

        if (next_lvl.empty()) break;
        list.insert(list.end(), next_lvl.begin(), next_lvl.end());
        curr = end;
    }

    return list;
}

void deleteFunnelTree(const vector<Funnel*> &list, const TriangleMesh& mesh, const indexType startIndex) {
    const size_t n = mesh.dictVertices[startIndex].size();
    for (size_t i = 0; i < n; i++) delete list[i];
    for (size_t i = n; i < list.size(); i+= 2) delete[] list[i];
}

TriangleMesh getMesh(const char *filename) {
    ifstream file(string("input/") + filename);
    size_t v, f, e;
    file >> v >> f >> e;

    vector<Point> points;
    points.reserve(v);
    for (size_t i = 0; i < v; i++) {
        double x, y, z;
        file >> x >> y >> z;
        points.emplace_back(x, y, z);
    }

    vector<Triangle> trianglesPointsIndexes;
    trianglesPointsIndexes.reserve(f);
    for (size_t i = 0; i < f; i++) {
        indexType a, b, c;
        short three;
        file >> three >> a >> b >> c;
        trianglesPointsIndexes.emplace_back(a, b, c);
    }

    return TriangleMesh(points, trianglesPointsIndexes);
}