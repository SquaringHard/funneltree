#define _USE_MATH_DEFINES   // M_PI
#include "ft.h"
#include <stdexcept>    // runtime_error
#include <algorithm>    // find, swap, min, max
#include <string>       // to_string
#include <cmath>        // acos, sqrt, INFINITY
#include <fstream>      // ifstream

TriangleMesh::TriangleMesh(const vector<Point> &ps, const vector<Triangle> &ts) : triangles(ts), points(ps), dictVertices(ps.size()) {
    if (ps.size() > MAX_INDEX) throw out_of_range("too many points");
    if (ts.size() > MAX_INDEX) throw out_of_range("too many faces");

    for (vector<Point>::iterator i = points.begin(); i != points.end(); i++) if (find(points.begin(), i, *i) != i)
        throw invalid_argument("point " + to_string(i - points.begin()) + " has duplicates");

    const indexType v = ps.size(), f = ts.size();
    dictVertices.reserve(v);
    dictEdges.reserve(v + f - 2);       // Euler's characteristic
    for (indexType i = 0; i < f; i++) { // add edges and triangle indexes to dictEdges
        for (const Edge &e : {Edge{ts[i][0], ts[i][1]}, Edge{ts[i][1], ts[i][2]}, Edge{ts[i][2], ts[i][0]}}) {
            const auto temp = dictEdges.try_emplace(e, array<indexType, 2>{i, MAX_INDEX});
            if (temp.second) continue;

            array<indexType, 2> &eFaces = temp.first->second;
            if (eFaces[1] != MAX_INDEX) throw logic_error("faces " + to_string(eFaces[0]) + ' ' + to_string(eFaces[1]) + ' ' + to_string(i) +
                                                          "occupying same edge");

            eFaces[1] = i;
        }

        for (const indexType j : ts[i]) dictVertices[j].push_back(i);
    }

    for (const pair<Edge, array<indexType, 2>> &temp : dictEdges) if (temp.second[1] == MAX_INDEX)
        { const Edge e = temp.first; throw logic_error("floating edge " + to_string(e[0]) + '-' + to_string(e[1])); }

    for (indexType i = 0; i < v; i++) if (dictVertices[i].empty()) throw logic_error("floating point " + to_string(i));
}

inline double angle(const double ab2, const double bc2, const double ca2) {
    const double cosABC = (ab2 + bc2 - ca2) / sqrt(ab2 * bc2) / 2;
    return cosABC >= 1 ? 0 : cosABC <= -1 ? M_PI : acos(cosABC);
}
inline double calPV2(const double pqv, const double pq2, const double qv2) { return pq2 + qv2 - sqrt(pq2 * qv2) * cos(pqv) * 2; }

double TriangleMesh::pistance2(const indexType a, const indexType b) const {
    const double abx = points[a].x - points[b].x, aby = points[a].y - points[b].y, abz = points[a].z - points[b].z;
    return abx * abx + aby * aby + abz * abz;
}

double TriangleMesh::pangle(const indexType a, const indexType b, const indexType c) const {
    const double abx = points[a].x - points[b].x, aby = points[a].y - points[b].y, abz = points[a].z - points[b].z,
                 cbx = points[c].x - points[b].x, cby = points[c].y - points[b].y, cbz = points[c].z - points[b].z;
    const double cosABC = (abx * cbx + aby * cby + abz * cbz) / sqrt((abx * abx + aby * aby + abz * abz) * (cbx * cbx + cby * cby + cbz * cbz));
    return cosABC >= 1 ? 0 : cosABC <= -1 ? M_PI : acos(cosABC);
}

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
        indexType p = psq[psq[0] == s], q = psq[(psq[2] == s) ? 1 : 2];

        double spq = mesh.pangle(s, p, q), psw = mesh.pangle(p, s, q);
        list[0] = Funnel(p, q, p, facesAt_s, mesh.pistance2(s, p), mesh.pistance2(p, q), spq, psw, psw, 0);
        swap(list[0].sequence[0], list[0].sequence[n - 1]);
        
        for (indexType i = 1; i < n; i++) {
            const array<indexType, 2> eFaces = mesh.dictEdges.at({s, q});
            psqIndex = eFaces[0] == psqIndex ? eFaces[1] : eFaces[0];
            psq = mesh.triangles[psqIndex];
            p = q; q = psq[0] + psq[1] + psq[2] - s - p;

            spq = mesh.pangle(s, p, q), psw = mesh.pangle(p, s, q);
            list[i] = Funnel(p, q, p, facesAt_s, mesh.pistance2(s, p), mesh.pistance2(p, q), spq, psw, psw, 0);
            vector<indexType> &sequence = list[i].sequence;
            swap(*find(sequence.begin(), sequence.end(), psqIndex), sequence[n - 1]);
        }
    }

    typedef unordered_map<Triangle, size_t, HashNComp> FunnelDict;
    FunnelDict twoChildrenFunnels;
    for (size_t start = 0, end = n;;) {
        const size_t curr_end = end;
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = start; i < curr_end; i++) {
            Funnel &f = list[i];
            if (f.removed) continue;

            find_v_suchthat_funnelhas2children:

            const array<indexType, 2> eFaces = mesh.dictEdges.at({f.x, f.q});
            const indexType nextFace = eFaces[eFaces[0] == f.sequence.back()];
            if (find(f.sequence.begin(), f.sequence.end(), nextFace) != f.sequence.end()) continue; // move to next i

            f.sequence.push_back(nextFace);
            indexType v;
            for (const indexType vNew : mesh.triangles[nextFace]) if (vNew != f.x && vNew != f.q) { v = vNew; break; }
            
            f.topright_angle += mesh.pangle(f.x, f.q, v);
            short int sign = 1;
            if (f.topright_angle >= M_PI) { f.topright_angle = M_PI * 2 - f.topright_angle; sign = -1; }
            const double vq2 = mesh.pistance2(v, f.q), pv2 = calPV2(f.topright_angle, f.pq2, vq2), spv = f.spq + angle(pv2, f.pq2, vq2) * sign;

            if (spv >= M_PI) { f.x = v; goto find_v_suchthat_funnelhas2children; }

            const double sv2 = calPV2(spv, f.sp2, pv2), psv = angle(f.sp2, sv2, pv2), pvq = angle(pv2, vq2, f.pq2), psw_new = min(f.psw, psv);
            f.topright_angle = max(mesh.pangle(f.x, v, f.q) - pvq * sign, 0.0);

            if (!(psv < f.psw)) { // bet ur temping to write psv >= psw
                f.q = v; f.pq2 = pv2; f.spq = spv; f.psq = psv; f.psw = psw_new;
                goto find_v_suchthat_funnelhas2children;
            }

            const double pvs = angle(pv2, sv2, f.sp2), vsq = f.psq - psv, vsw = f.psw - psv, svq = pvq - pvs;
            f.pvs = pvs;
            #pragma omp atomic capture
            { f.childrenIndex = end; end += 2; }
            list[f.childrenIndex] = Funnel(f.p, v, f.x, f.sequence, f.sp2, pv2, spv, psv, psw_new, f.topright_angle);
            list[f.childrenIndex + 1] = Funnel(v, f.q, v, f.sequence, sv2, vq2, svq, vsq, vsw, 0);

            pair<FunnelDict::iterator, bool> temp;
            #pragma omp critical
            temp = twoChildrenFunnels.try_emplace({f.p, v, f.q}, i);
            if (temp.second) continue;

            const Funnel &oldFunnel = list[temp.first->second];
            const size_t oldChildrenIndex = oldFunnel.childrenIndex;
            const double oldsv2 = list[oldChildrenIndex + 1].sp2;
            const bool pvs_is_larger_than_pvs2 = pvs > oldFunnel.pvs;

            if (oldsv2 > sv2) { flag(list, oldChildrenIndex + pvs_is_larger_than_pvs2); temp.first->second = i; }
            else if (sv2 > oldsv2) list[f.childrenIndex + !pvs_is_larger_than_pvs2].removed = true;
            else { list[f.childrenIndex + !pvs_is_larger_than_pvs2].removed = true; flag(list, oldChildrenIndex + pvs_is_larger_than_pvs2); }
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