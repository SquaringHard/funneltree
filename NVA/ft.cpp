#define _USE_MATH_DEFINES   // M_PI
#include "ft.h"
#include <stdexcept>    // runtime_error
#include <iterator>     // distance
#include <string>
#include <limits>       // numeric_limits
#include <cmath>        // acos, sqrt, INFINITY
#include <algorithm>    // find, swap
#include <omp.h>

constexpr const indexType MAX_INDEX = numeric_limits<indexType>::max();

TriangleMesh::TriangleMesh(const vector<Point> &ps, const vector<Triangle> &ts) : triangles(ts), points(ps), dictVertices(ps.size()) {
    const size_t v = ps.size(), f = ts.size();
    if (v > MAX_INDEX) throw runtime_error("too many points");
    if (f > MAX_INDEX) throw runtime_error("too many faces");
    dictEdges.reserve(v + f - 2);   // Euler's characteristic

    for (auto i = points.begin(); i != points.end(); i++)
        if (find(points.begin(), i, *i) != i) throw runtime_error("point " + to_string(distance(i, points.begin())) + " has duplicates");

    for (indexType i = 0; i < f; i++) {
        for (const Edge &e : {Edge(ts[i].a, ts[i].b), Edge(ts[i].b, ts[i].c), Edge(ts[i].c, ts[i].a)}) {    // add edges and triangle indexes to dictEdges
            const pair<DictEdgeType::iterator, bool> temp = dictEdges.try_emplace(e, array<indexType, 2>{i, MAX_INDEX});
            if (temp.second) continue;
            array<indexType, 2> &eFaces = temp.first->second;
            if (eFaces[1] != MAX_INDEX) throw runtime_error("faces " + to_string(eFaces[0]) + ' ' + to_string(eFaces[1]) + ' ' + to_string(i) + "occupying same edge");
            eFaces[1] = i;
        }

        for (const indexType j : {ts[i].a, ts[i].b, ts[i].c}) dictVertices[j].push_back(i);   // add triangle indexes to dictVertices
    }

    for (const pair<Edge, array<indexType, 2>> &temp : dictEdges) if (temp.second[1] == MAX_INDEX) {
            const Edge e = temp.first;
            throw runtime_error("floating edge " + to_string(e.a) + '-' + to_string(e.b));
    }

    for (indexType i = 0; i < v; i++) if (dictVertices[i].empty()) throw runtime_error("floating point " + to_string(i));
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

vector<double> tol(omp_get_max_threads());

vector<Funnel*> FunnelTree(const TriangleMesh& mesh, const indexType s) {
    const vector<indexType> &facesAt_s = mesh.dictVertices[s];
    vector<Funnel*> list;
    for (indexType i = 0; i < facesAt_s.size(); i++) {  // if parallelizing this, read deleteFunnelTree() first
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

    #pragma omp declare reduction(merge: vector<Funnel*>: omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    typedef unordered_map<Triangle, const Funnel*, HashNComp, HashNComp> FunnelDict;
    FunnelDict twoChildrenFunnels;
    for (size_t curr_lvl_index = 0, next_lvl_index = list.size();;) {
        vector<Funnel*> next_lvl;
        #pragma omp parallel for reduction(merge: next_lvl) schedule(dynamic)
        for (size_t i = curr_lvl_index; i < next_lvl_index; i++) {
            const double start = omp_get_wtime();
            Funnel *const funnel = list[i];
            if (funnel->removed) 
                {tol[omp_get_thread_num()] += omp_get_wtime() - start; continue;}

            const indexType p = funnel->p;
            indexType q = funnel->q, x = funnel->x, v;
            const double sp = funnel->sp;
            double topright_angle = funnel->topright_angle, pq = funnel->pq, pv, vq, spv, psv, pvq, psw, sv;
            vector<indexType> &sequence = funnel->sequence;

            while (true) {
                spv = INFINITY;
                short int sign;
                while (true) {
                    const array<indexType, 2> eFaces = mesh.dictEdges.at(Edge(x, q));
                    const indexType nextFace = eFaces[0] == sequence.back() ? eFaces[1] : eFaces[0];
                    if (find(sequence.begin(), sequence.end(), nextFace) != sequence.end()) break;

                    sequence.push_back(nextFace);
                    const Triangle temp2 = mesh.triangles[nextFace];
                    for (const indexType vNew : {temp2.a, temp2.b, temp2.c}) if (vNew != x && vNew != q) { v = vNew; break; }
                    
                    topright_angle += mesh.pangle(x, q, v);
                    sign = 1;
                    if (topright_angle >= M_PI) { topright_angle = M_PI * 2 - topright_angle; sign = -1; }

                    vq = mesh.pistance(v, q);
                    pv = calPV(topright_angle, pq, vq);
                    spv = funnel->spq + angle(pv, pq, vq) * sign;

                    if (spv >= M_PI) x = v;
                    else break;
                }

                if (spv >= M_PI) break;
                sv = calPV(spv, sp, pv);
                psv = angle(sp, sv, pv);
                pvq = angle(pv, vq, pq);
                psw = min(funnel->psw, psv);
                topright_angle = max(mesh.pangle(x, v, q) - pvq * sign, 0.0);

                if (psv < funnel->psw) break;
                q = v;
                pq = pv;
                funnel->spq = spv;
                funnel->psq = psv;
                funnel->psw = psw;
            }

            if (spv >= M_PI) {tol[omp_get_thread_num()] += omp_get_wtime() - start; continue;}
            
            const double pvs = angle(pv, sv, sp), vsq = funnel->psq - psv, vsw = funnel->psw - psv;
            funnel->pvs = pvs;
            funnel->children = new Funnel[2]{Funnel(p, v, x, sequence, sp, pv, spv, psv, psw, topright_angle),
                                             Funnel(v, q, v, sequence, sv, vq, pvq - pvs, vsq, vsw)};
            Funnel *const child0 = funnel->children;
            next_lvl.push_back(child0);
            next_lvl.push_back(child0 + 1);

            pair<FunnelDict::iterator, bool> temp;
            #pragma omp critical
            temp = twoChildrenFunnels.try_emplace({p, v, q}, funnel);
            if (temp.second) {tol[omp_get_thread_num()] += omp_get_wtime() - start; continue;}

            /////////////////////////// CLIP OFF FUNNELS ////////////////////////////

            const Funnel *&oldFunnel = temp.first->second;
            Funnel *const oldChild0 = oldFunnel->children;
            const double sv2 = (oldChild0 + 1)->sp, pvs2 = oldFunnel->pvs;

            if (sv2 > sv) {
                if (pvs > pvs2) (oldChild0 + 1)->remove();
                else oldChild0->remove();
                oldFunnel = funnel; // switch marked funnel with current funnel, which has shorter sv length and therefore may help remove more funnels
            } else if (sv > sv2) {
                if (pvs > pvs2) child0->removed = true;
                else (child0 + 1)->removed = true;
            } else if (pvs > pvs2) {
                child0->removed = true;
                (oldChild0 + 1)->remove();
            } else {
                (child0 + 1)->removed = true;
                oldChild0->remove();
            }

            tol[omp_get_thread_num()] += omp_get_wtime() - start;
        }

        if (next_lvl.empty()) break;
        list.insert(list.end(), next_lvl.begin(), next_lvl.end());
        curr_lvl_index = next_lvl_index;
        next_lvl_index = list.size();
    }

    return list;
}

void deleteFunnelTree(const vector<Funnel*> &list, const TriangleMesh& mesh, const indexType startIndex) {
    const size_t n = mesh.dictVertices[startIndex].size();
    for (size_t i = 0; i < n; i++) delete list[i];
    for (size_t i = n; i < list.size(); i+= 2) delete[] list[i];
}