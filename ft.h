#ifndef FT_H
#define FT_H
using namespace std;
#include <vector>
#include <unordered_map>
#include <array>
#define THREAD_TIMING


typedef int indexType;
const size_t MAX_INDEX = 1e8;   // max number of vertices for HashNComp()(Triangle) to work properly

#ifdef THREAD_TIMING
    #include <chrono>
    extern vector<chrono::nanoseconds> threadRuntime, threadIdleTime;
#endif

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
    size_t operator()(const Triangle &t) const {
        const size_t seed = 0x9e3779b * (1 + t.a) + t.b;
        return ((seed << 27) + 0x517CC1B7 * (t.c)) ^ (seed >> 37);
    }
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
TriangleMesh getMesh(const char *filename, const indexType startIndex);


#endif