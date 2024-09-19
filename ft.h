#ifndef FT_H
#define FT_H
using namespace std;
#include <array>
#include <vector>
#include <unordered_map>


typedef int indexType;
const size_t MAX_INDEX = 1e8;   // max number of vertices for HashNComp()(Triangle) to work properly

struct Point {
    double x, y, z;
    bool operator==(const Point &p) const { return x == p.x && y == p.y && z == p.z; }  // for find() in TriangleMesh constructor
};

typedef array<indexType, 2> Edge;
typedef array<indexType, 3> Triangle;

struct HashNComp {
    size_t operator()(const Edge &e) const { return e[0] ^ e[1]; }
    bool operator()(const Edge &a, const Edge &b) const { return a[0] == b[0] && a[1] == b[1] || a[0] == b[1] && a[1] == b[0]; }
    size_t operator()(const Triangle &t) const {
        const size_t seed = 0x9e3779b * (1 + t[0]) + t[1];
        return ((seed << 27) + 0x517CC1B7 * (t[2])) ^ (seed >> 37);
    }
};

struct TriangleMesh {
    const vector<Triangle> triangles;
    const vector<Point> points;
    vector<vector<indexType>> dictVertices;
    unordered_map<Edge, array<indexType, 2>, HashNComp, HashNComp> dictEdges;
    TriangleMesh(const vector<Point> &points, const vector<Triangle> &triangles);
    double pistance2(const indexType a, const indexType b) const;
    double pangle(const indexType a, const indexType b, const indexType c) const;
};

struct Funnel {
    size_t childrenIndex;
    vector<indexType> sequence;
    double sp2, pq2, spq, psq, psw, topright_angle, pvs;
    indexType p, q, x;
    bool removed;
    Funnel() = default;
    Funnel(const indexType p, const indexType q, const indexType x, const vector<indexType> &sequence, const double sp2, const double pq2,
           const double spq, const double psq, const double psw, const double topright_angle)
    : childrenIndex(0), sequence(sequence), p(p), q(q), x(x), sp2(sp2), pq2(pq2), spq(spq), psq(psq), psw(psw), topright_angle(topright_angle), removed(false) {}
};

vector<Funnel> FunnelTree(const TriangleMesh& mesh, const indexType startIndex);
TriangleMesh getMesh(const char *filename);


#endif