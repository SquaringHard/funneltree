#ifndef FT_H
#define FT_H
using namespace std;
#include <functional>   // hash function
#include <array>
#include <vector>
#include <utility>      // pair

typedef int indexType;

struct Point {
    double x, y, z;
    Point(const double x, const double y, const double z) : x(x), y(y), z(z) {}
    bool operator==(const Point &p) const { return x == p.x && y == p.y && z == p.z; }  // for find() in TriangleMesh constructor
};

struct Edge {
    indexType a, b;
    Edge(const indexType a, const indexType b) : a(a), b(b) {}
};

struct Triangle {
    indexType a, b, c;
    Triangle(const indexType a, const indexType b, const indexType c) : a(a), b(b), c(c) {}
};

struct Funnel {
    Funnel *children;
    vector<indexType> sequence;
    double sp, pq, spq, psq, psw, topright_angle, pvs;
    indexType p, q, x;
    bool removed;
    Funnel(const indexType p, const indexType q, const indexType x, const vector<indexType> &sequence, const double sp, const double pq,
           const double spq, const double psq, const double psw, const double topright_angle = 0)
    : children(nullptr), sequence(sequence), p(p), q(q), x(x), sp(sp), pq(pq), spq(spq), psq(psq), psw(psw), topright_angle(topright_angle), removed(false) {}
    void remove();
};

struct HashNComp {
    size_t operator()(const Edge &e) const { return e.a ^ e.b; }
    size_t operator()(const Triangle &t) const { return ((t.a ^ (t.b << 1)) >> 1) ^ (t.c << 1); }
    bool operator()(const Edge &a, const Edge &b) const { return a.a == b.a && a.b == b.b || a.a == b.b && a.b == b.a; }
    bool operator()(const Triangle &a, const Triangle &b) const { return a.a == b.a && a.b == b.b && a.c == b.c; }
};

struct TriangleMesh {
    const vector<Triangle> triangles;
    const vector<Point> points;
    vector<vector<indexType>> dictVertices;
    typedef unordered_map<Edge, array<indexType, 2>, HashNComp, HashNComp> DictEdgeType;    // an edge is on exactly 2 faces
    DictEdgeType dictEdges;
    TriangleMesh(const vector<Point> &points, const vector<Triangle> &triangles);
    double pistance(const indexType a, const indexType b) const;
    double pangle(const indexType a, const indexType b, const indexType c) const;
};

vector<Funnel*> FunnelTree(const TriangleMesh& mesh, const indexType startIndex);
void deleteFunnelTree(const vector<Funnel*> &list, const TriangleMesh& mesh, const indexType startIndex);


#endif