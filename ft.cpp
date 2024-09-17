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
        for (const Edge &e : {Edge{ts[i].a, ts[i].b}, Edge{ts[i].b, ts[i].c}, Edge{ts[i].c, ts[i].a}}) {
            const auto temp = dictEdges.try_emplace(e, array<indexType, 2>{i, MAX_INDEX});
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

void flag(vector<Funnel> &list, const size_t i) {
    Funnel &f = list[i];
    if (f.childrenIndex == 0) { f.removed = true; return; }
    flag(list, f.childrenIndex);
    flag(list, f.childrenIndex + 1);
}

vector<Funnel> FunnelTree(const TriangleMesh& mesh, const indexType s) {
    const vector<indexType> &facesAt_s = mesh.dictVertices[s];
    const indexType n = facesAt_s.size();
    vector<Funnel> list(n * 3); // n + n * 2

    {
        indexType psqIndex = facesAt_s[0];
        Triangle psq = mesh.triangles[psqIndex];
        indexType p = psq.a == s ? psq.b : psq.a, q = psq.a + psq.b + psq.c - s - p;

        double spq = mesh.pangle(s, p, q), psw = mesh.pangle(p, s, q);
        list[0] = Funnel(p, q, p, facesAt_s, mesh.pistance(s, p), mesh.pistance(p, q), spq, psw, psw, 0);
        swap(list[0].sequence[0], list[0].sequence[n - 1]);
        
        for (indexType i = 1; i < n; i++) {
            const array<indexType, 2> eFaces = mesh.dictEdges.at({s, q});
            psqIndex = eFaces[0] == psqIndex ? eFaces[1] : eFaces[0];
            psq = mesh.triangles[psqIndex];
            p = q; q = psq.a + psq.b + psq.c - s - p;

            spq = mesh.pangle(s, p, q), psw = mesh.pangle(p, s, q);
            list[i] = Funnel(p, q, p, facesAt_s, mesh.pistance(s, p), mesh.pistance(p, q), spq, psw, psw, 0);
            vector<indexType> &sequence = list[i].sequence;
            swap(*find(sequence.begin(), sequence.end(), psqIndex), sequence[n - 1]);
        }
    }

    typedef unordered_map<Triangle, size_t, HashNComp, HashNComp> FunnelDict;
    FunnelDict twoChildrenFunnels;
    for (size_t start = 0, end = n;;) {
        const size_t curr_end = end;
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = start; i < curr_end; i++) {
            Funnel &f = list[i];
            if (f.removed) continue;
            indexType p = f.p; double sp = f.sp;

            find_v_suchthat_funnelhas2children:

            const array<indexType, 2> eFaces = mesh.dictEdges.at({f.x, f.q});
            const indexType nextFace = eFaces[0] == f.sequence.back() ? eFaces[1] : eFaces[0];
            if (find(f.sequence.begin(), f.sequence.end(), nextFace) != f.sequence.end()) continue; // move to next i

            f.sequence.push_back(nextFace);
            indexType v;
            const Triangle temp2 = mesh.triangles[nextFace];
            for (const indexType vNew : {temp2.a, temp2.b, temp2.c}) if (vNew != f.x && vNew != f.q) { v = vNew; break; }
            
            f.topright_angle += mesh.pangle(f.x, f.q, v);
            short int sign = 1;
            if (f.topright_angle >= M_PI) { f.topright_angle = M_PI * 2 - f.topright_angle; sign = -1; }
            const double vq = mesh.pistance(v, f.q), pv = calPV(f.topright_angle, f.pq, vq), spv = f.spq + angle(pv, f.pq, vq) * sign;

            if (spv >= M_PI) { f.x = v; goto find_v_suchthat_funnelhas2children; }

            const double sv = calPV(spv, sp, pv), psv = angle(sp, sv, pv), pvq = angle(pv, vq, f.pq), psw_new = min(f.psw, psv);
            f.topright_angle = max(mesh.pangle(f.x, v, f.q) - pvq * sign, 0.0);

            if (!(psv < f.psw)) { // bet ur temping to write psv >= psw
                make_only_Fpv:
                    f.q = v; f.pq = pv; f.spq = spv; f.psq = psv; f.psw = psw_new;
                    goto find_v_suchthat_funnelhas2children;
            }

            const double pvs = angle(pv, sv, sp), vsq = f.psq - psv, vsw = f.psw - psv, svq = pvq - pvs;
            FunnelDict::const_iterator temp;
            #pragma omp critical (access_pvq)
            temp = twoChildrenFunnels.find({p, v, f.q});
            if (temp != twoChildrenFunnels.end()) { // clip_off_funnel
                const Funnel &oldFunnel = list[temp->second];
                const size_t oldChildrenIndex = oldFunnel.childrenIndex;
                const double sv2 = list[oldChildrenIndex + 1].sp;
                const bool pvs_is_larger_than_pvs2 = pvs > oldFunnel.pvs;   // i feel like it

                if (sv2 > sv) flag(list, oldChildrenIndex + pvs_is_larger_than_pvs2);
                else {
                    if (sv > sv2) { if (!pvs_is_larger_than_pvs2) goto make_only_Fpv; }
                    else if (pvs_is_larger_than_pvs2) flag(list, oldChildrenIndex + 1);
                    else { flag(list, oldChildrenIndex); goto make_only_Fpv; }

                    // make only Fvq
                    p = v; f.x = v; sp = sv; f.pq = vq; f.spq = svq; f.psq = vsq; f.psw = vsw; f.topright_angle = 0;
                    goto find_v_suchthat_funnelhas2children;
                }
            }

            // make Fpv and Fvq:
            f.pvs = pvs;
            #pragma omp atomic capture
            { f.childrenIndex = end; end += 2; }
            list[f.childrenIndex] = Funnel(p, v, f.x, f.sequence, sp, pv, spv, psv, psw_new, f.topright_angle);
            list[f.childrenIndex + 1] = Funnel(v, f.q, v, f.sequence, sv, vq, svq, vsq, vsw, 0);
            #pragma omp critical (access_pvq)
            twoChildrenFunnels[{p, v, f.q}] = i;
        }

        if (end == curr_end) { list.resize(end); break; }
        start = curr_end;
        list.resize((end - start) * 3 + start);
    }

    return list;
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
        points.push_back({x, y, z});
    }

    vector<Triangle> trianglesPointsIndexes;
    trianglesPointsIndexes.reserve(f);
    for (size_t i = 0; i < f; i++) {
        indexType a, b, c;
        short three;
        file >> three >> a >> b >> c;
        trianglesPointsIndexes.push_back({a, b, c});
    }

    return TriangleMesh(points, trianglesPointsIndexes);
}