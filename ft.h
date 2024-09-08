#ifndef FT_H
#define FT_H
using namespace std;
#include <functional>   // hash function
#include <array>
#include <vector>

typedef int indexType;

struct Point {  // use pointers because they seems to be faster
    double x, y, z;
    indexType index;
    Point(const double x, const double y, const double z, const indexType index) : x(x), y(y), z(z), index(index) {}
    bool operator==(const Point &p) const { return x == p.x && y == p.y && z == p.z; }
};

struct Edge {
    const Point *a, *b;
    Edge(const Point *const a, const Point *const b) : a(a), b(b) {}
    bool operator==(const Edge &e) const { return *a == *e.a && *b == *e.b || *a == *e.b && *b == *e.a; }
};

struct Triangle {
    const Point *a, *b, *c;
    Triangle(const Point *const a, const Point *const b, const Point *const c) : a(a), b(b), c(c) {}
};

struct Hasher {
    size_t operator()(const Point &p) const { return ((hash<double>()(p.x) ^ (hash<double>()(p.y) << 1)) >> 1) ^ (hash<double>()(p.z) << 1); }
    size_t operator()(const Edge &e) const { return Hasher()(*e.a) ^ Hasher()(*e.b); }
    size_t operator()(const array<const Point*, 3> &t) const { return ((Hasher()(*t[0]) ^ (Hasher()(*t[1]) << 1)) >> 1) ^ (Hasher()(*t[2]) << 1); }
};

typedef unordered_map<Point, vector<indexType>, Hasher> PointDict;
typedef unordered_map<Edge, array<indexType, 2>, Hasher> EdgeDict;    // an edge is on exactly 2 faces

struct TriangleMesh {
    vector<Triangle> triangles;
    PointDict dictVertices;
    EdgeDict dictEdges;
    TriangleMesh(const vector<Point> &points, const vector<array<indexType, 3>> &trianglesPointsIndexes);
};

struct Funnel {
    Funnel *childPV, *childVQ;
    vector<indexType> sequence;
    const Point *p, *q, *x;
    double sp, pq, spq, psq, psw, topright_angle;
    bool removed;
    Funnel(const Point *const p, const Point *const q, const Point *const x, const vector<indexType> &sequence, const double sp, const double pq,
           const double spq, const double psq, const double psw, const double topright_angle = 0)
    : childPV(nullptr), childVQ(nullptr), sequence(sequence), p(p), q(q), x(x), sp(sp), pq(pq), spq(spq), psq(psq), psw(psw), topright_angle(topright_angle), removed(false) {}
    void remove();
};

vector<Funnel*> FunnelTree(const TriangleMesh& mesh, const Point &s);
inline void deleteFunnelTree(const vector<Funnel*> &list) { for (Funnel *f : list) delete f; }


#endif