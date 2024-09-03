// run all files: ./main
// run some files: ./main cube1.geom sphere4.geom
// to measure runtime, replace run function in main with time function

#include "ft.h"
#include <iostream>     // cin, cout
#include <string>
#include <fstream>      // ifstream
#include <chrono>
#include <cmath>        // INFINITY
#include <iomanip>      // setprecision
#include <algorithm>    // equal
#include <numeric>      // reduce
#define COMPARE_LENGTH_MODE

bool compareFiles(const char *filename) { // https://stackoverflow.com/a/37575457
    ifstream f1(string("expectedLength/") + filename, ifstream::binary|ifstream::ate);
    ifstream f2(string("outputLength/") + filename, ifstream::binary|ifstream::ate);

    if (f1.fail() || f2.fail()) return false;     // file problem
    if (f1.tellg() != f2.tellg()) return false;   // size mismatch

    // seek back to beginning and use equal to compare contents
    f1.seekg(0, ifstream::beg);
    f2.seekg(0, ifstream::beg);
    return equal(istreambuf_iterator<char>(f1.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(f2.rdbuf()));
}

void run(const char *filename) {
    cout << "File \"" << filename << "\": ";
    ifstream file(string("input/") + filename);

    size_t v, f, e;
    file >> v >> f >> e;

    vector<Point> points;
    points.reserve(v);
    for (size_t i = 0; i < v; i++) {
        double x, y, z;
        file >> x >> y >> z;
        points.emplace_back(x, y, z, i);
    }

    vector<array<indexType, 3>> trianglesPointsIndexes;
    trianglesPointsIndexes.reserve(f);
    for (size_t i = 0; i < f; i++) {
        indexType a, b, c, three;
        file >> three >> a >> b >> c;
        trianglesPointsIndexes.push_back({a, b, c});
    }

    TriangleMesh mesh(points, trianglesPointsIndexes);
    // cout << "\nThere are " << mesh.dictVertices.size() << " vertices, " << mesh.triangles.size() << " faces and " << mesh.dictEdges.size() << " edges.\n";
    
    const size_t i = 0;
    const auto start = chrono::high_resolution_clock::now();
    const vector<Funnel*> list = FunnelTree(points[i], mesh);
    const auto end = chrono::high_resolution_clock::now();
    const auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Funnel tree initialized with " << list.size() << " nodes in " << duration.count() << " milliseconds.\n";

    #ifdef COMPARE_LENGTH_MODE
        vector<double> shortestLength(v, INFINITY);
        shortestLength[i] = 0;
        for (const Funnel *const f : list)
            if (shortestLength[f->p->index] > f->sp) shortestLength[f->p->index] = f->sp;
        ofstream out(string("outputLength/") + filename);
        out << fixed << setprecision(8);
        for (const double i : shortestLength) out << i << '\n';
        out.flush();
        if (!compareFiles(filename)) cout << "---------- NOT PASSED ----------\n";
    #endif

    for (const Funnel *funnel : list) delete funnel;
}

void time(const char *filename, const size_t n = 100) {
    cout << "File \"" << filename << "\" ran " << n << " times. Avg: ";
    ifstream file(string("input/") + filename);

    size_t v, f, e;
    file >> v >> f >> e;

    vector<Point> points;
    points.reserve(v);
    for (size_t i = 0; i < v; i++) {
        double x, y, z;
        file >> x >> y >> z;
        points.emplace_back(x, y, z, i);
    }

    vector<array<indexType, 3>> trianglesPointsIndexes;
    trianglesPointsIndexes.reserve(f);
    for (size_t i = 0; i < f; i++) {
        indexType a, b, c, three;
        file >> three >> a >> b >> c;
        trianglesPointsIndexes.push_back({a, b, c});
    }

    TriangleMesh mesh(points, trianglesPointsIndexes);

    vector<chrono::nanoseconds> duration;
    for (size_t i = 0; i < n; i++) {
        const auto start = chrono::high_resolution_clock::now();
        const vector<Funnel*> list = FunnelTree(points[0], mesh);
        const auto end = chrono::high_resolution_clock::now();
        for (const Funnel *funnel : list) delete funnel;
        duration.emplace_back(end - start);
    }

    const auto avg = reduce(duration.begin(), duration.end()) / n;
    chrono::nanoseconds error(0);
    
    // rand error
    for (const auto i : duration) error += i < avg ? avg - i : i - avg;
    error /= n;

    error++;    // sys error
    cout << chrono::duration_cast<chrono::microseconds>(avg).count() << " +/- "
         << chrono::duration_cast<chrono::microseconds>(error).count() << " microseconds.\n";
}

int main(int argc, const char *argv[]) {
    const char *allfiles[] = {"main", "cube1.geom", "cube2.geom", "cube3.geom", "cube4.geom",
                              "sphere1.geom", "sphere2.geom", "sphere3.geom", "sphere4.geom",
                              "spiral1.geom", "spiral2.geom", "J17.geom"},
               **files;
    if (argc > 1) files = argv; else { files = allfiles; argc = sizeof(allfiles) / sizeof(*allfiles); }
    for (int i = 1; i < argc; i++) {
        // run(files[i]);
        time(files[i]);
    }
}