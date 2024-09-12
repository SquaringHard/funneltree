#define _USE_MATH_DEFINES   // M_PI
#include "ft.h"
#include <stdexcept>    // runtime_error
#include <string>
#include <limits>       // numeric_limits
#include <cmath>        // acos, sqrt, INFINITY
#include <algorithm>    // find

constexpr const indexType MAX_INDEX = numeric_limits<indexType>::max();

TriangleMesh::TriangleMesh(const vector<Point> &points, const vector<array<indexType, 3>> &trianglesPointsIndexes) {
    const size_t v = points.size(), f = trianglesPointsIndexes.size();
    if (v > MAX_INDEX) throw runtime_error("too many points");
    if (f > MAX_INDEX) throw runtime_error("too many faces");
    dictVertices.reserve(v);
    dictEdges.reserve(v + f - 2);   // Euler's characteristic

    vector<const Point*> pointsPointers;
    for (const Point &p : points) { // add points
        const pair<PointDict::iterator, bool> temp = dictVertices.try_emplace(p, vector<indexType>());
        if (!temp.second) throw runtime_error("point " + to_string(p.index) + " has duplicates");
        pointsPointers.push_back(&temp.first->first);
    }

    for (indexType i = 0; i < f; i++) {
        const Point *const a = pointsPointers[trianglesPointsIndexes[i][0]], *const b = pointsPointers[trianglesPointsIndexes[i][1]], *const c = pointsPointers[trianglesPointsIndexes[i][2]];

        for (const Edge &e : {Edge(a, b), Edge(b, c), Edge(c, a)}) {            // add edges and triangle indexes to dictEdges
            const pair<EdgeDict::iterator, bool> temp = dictEdges.try_emplace(e, array<indexType, 2>{i, MAX_INDEX});
            if (temp.second) continue;
            array<indexType, 2> &eFaces = temp.first->second;
            if (eFaces[1] != MAX_INDEX) throw runtime_error("faces " + to_string(eFaces[0]) + ' ' + to_string(eFaces[1]) + ' ' + to_string(i) + "occupying same edge");
            eFaces[1] = i;
        }

        for (const Point *const p : {a, b, c}) dictVertices[*p].push_back(i);   // add triangle indexes to dictVertices
        triangles.emplace_back(a, b, c);                                        // add triangles
    }

    for (const pair<Edge, array<indexType, 2>> &temp : dictEdges) if (temp.second[1] == MAX_INDEX) {
            const Edge e = temp.first;
            throw runtime_error("floating edge " + to_string(e.a->index) + '-' + to_string(e.b->index));
    }

    for (const pair<Point, vector<indexType>> &temp : dictVertices) if (temp.second.empty()) throw runtime_error("floating point " + to_string(temp.first.index));
}

inline double distance(const Point *const a, const Point *const b) {
    const double abx = a->x - b->x, aby = a->y - b->y, abz = a->z - b->z;
    return sqrt(abx * abx + aby * aby + abz * abz);
}

inline double angle(const double ab, const double bc, const double ca) { return acos((ab * ab + bc * bc - ca * ca) / (ab * bc * 2)); }
inline double angle(const Point *const a, const Point *const b, const Point *const c) { return angle(distance(a, b), distance(b, c), distance(c, a)); }
inline double calPV(const double pqv, const double pq, const double qv) { return sqrt(pq * pq + qv * qv - pq * qv * cos(pqv) * 2); }

void Funnel::remove() {
    if (childPV == nullptr) { removed = true; return; }
    childPV->remove();
    childVQ->remove();
}

struct FunnelHashCompare {
    size_t operator()(const array<const Point*, 3> &t) const { return ((Hasher()(*t[0]) ^ (Hasher()(*t[1]) << 1)) >> 1) ^ (Hasher()(*t[2]) << 1); }
    bool operator()(const array<const Point*, 3> &a, const array<const Point*, 3> &b) const { return *a[0] == *b[0] && *a[1] == *b[1] && *a[2] == *b[2]; }
};

typedef unordered_map<array<const Point*, 3>, const Funnel*, FunnelHashCompare, FunnelHashCompare> FunnelDict;

vector<Funnel*> FunnelTree(const TriangleMesh& mesh, const Point &s) {
    const PointDict::const_iterator dictVerticesAt_s = mesh.dictVertices.find(s);
    if (dictVerticesAt_s == mesh.dictVertices.end()) return {};
    const vector<indexType> &facesAt_s = dictVerticesAt_s->second;

    vector<Funnel*> list;
    for (indexType i = 0; i < facesAt_s.size(); i++) {
        const Triangle pqv = mesh.triangles[facesAt_s[i]];
        const Point *p, *q;
        if (s == *pqv.a) {
            q = pqv.b;
            p = pqv.c;
        } else if (s == *pqv.b) {
            q = pqv.c;
            p = pqv.a;
        } else {
            q = pqv.a;
            p = pqv.b;
        }

        const double spq = angle(&s, p, q), psw = angle(p, &s, q);

        // insert all faces containing s instead of just 1 face to sequence so that the funnels never reach these faces again
        // because then the funnels don't have to reach the vertices on these faces
        vector<indexType> newSequence = facesAt_s;
        swap(newSequence[i], newSequence.back());
        list.push_back(new Funnel(p, q, p, newSequence, distance(&s, p), distance(p, q), spq, psw, psw));
    }

    FunnelDict twoChildrenFunnels;
    for (size_t i = 0; i < list.size(); i++) {
        Funnel *const funnel = list[i];
        if (funnel->removed) continue;

        const Point *const p = funnel->p, *q = funnel->q, *x = funnel->x, *v;
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
                for (const Point *const vNew : {temp2.a, temp2.b, temp2.c}) if (!(*vNew == *x || *vNew == *q)) { v = vNew; break; }
                
                topright_angle += angle(x, q, v);
                sign = 1;
                if (topright_angle >= M_PI) { topright_angle = M_PI * 2 - topright_angle; sign = -1; }

                vq = distance(v, q);
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
            topright_angle = max(angle(x, v, q) - pvq * sign, 0.0);

            if (psv < funnel->psw) break;
            q = v;
            pq = pv;
            funnel->spq = spv;
            funnel->psq = psv;
            funnel->psw = psw;
        }

        if (spv >= M_PI) continue;
        
        const double pvs = angle(pv, sv, sp), vsq = funnel->psq - psv, vsw = funnel->psw - psv;
        Funnel *const childPV = new Funnel(p, v, x, sequence, sp, pv, spv, psv, psw, topright_angle),
               *const childVQ = new Funnel(v, q, v, sequence, sv, vq, pvq - pvs, vsq, vsw);
        funnel->childPV = childPV;
        funnel->childVQ = childVQ;
        list.push_back(childPV);
        list.push_back(childVQ);

        const auto temp = twoChildrenFunnels.try_emplace({p, v, q}, funnel);
        if (temp.second) continue;
        
        /////////////////////////// CLIP OFF FUNNELS ////////////////////////////

        const Funnel *&oldFunnel = temp.first->second;
        Funnel *const oldChildPV = oldFunnel->childPV, *const oldChildVQ = oldFunnel->childVQ;
        const double sv2 = oldChildVQ->sp, pvs2 = asin(oldChildPV->sp * sin(oldChildPV->spq) / sv2);

        if (sv2 > sv) {
            if (pvs > pvs2) oldChildVQ->remove();
            else oldChildPV->remove();
            oldFunnel = funnel; // switch marked funnel with current funnel, which has shorter sv length and therefore may help remove more funnels
        } else if (sv > sv2) {
            if (pvs > pvs2) childPV->removed = true;
            else childVQ->removed = true;
        } else if (pvs > pvs2) {
            childPV->removed = true;
            oldChildVQ->remove();
        } else {
            childVQ->removed = true;
            oldChildPV->remove();
        }
    }

    return list;
}