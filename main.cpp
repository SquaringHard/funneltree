// run all files: ./main
// run some files: ./main cube1.geom sphere4.geom
// to measure runtime, replace run function in main with time function

#include "ft.h"
#include <utility>      // pair
#include <string>
#include <fstream>      // ifstream, ofstream
#include <cmath>        // INFINITY, fabs
#include <filesystem>
#include <iomanip>      // fixed, setprecision
#include <limits>       // numeric_limits
#include <iostream>     // cin, cout
#include <chrono>
#include <algorithm>    // equal
#include <numeric>      // reduce

pair<TriangleMesh, Point> getMeshAndStartPoint(const char *filename, const size_t startIndex) {
    ifstream file(string("input/") + filename);
    size_t v, f, e;
    file >> v >> f >> e;

    vector<Point> points;
    for (size_t i = 0; i < v; i++) {
        double x, y, z;
        file >> x >> y >> z;
        points.emplace_back(x, y, z, i);
    }

    vector<array<indexType, 3>> trianglesPointsIndexes;
    for (size_t i = 0; i < f; i++) {
        indexType a, b, c;
        short three;
        file >> three >> a >> b >> c;
        trianglesPointsIndexes.push_back({a, b, c});
    }

    return {TriangleMesh(points, trianglesPointsIndexes), points.at(startIndex)};
}

bool compareLength(const char *filename, const vector<Funnel*> &list, const size_t startIndex, const size_t n = 0) {
    const string realFilename = string(filename) + "_s=" + to_string(startIndex);
    ifstream file(string("expected/") + realFilename + ".txt");
    if (!file.is_open()) throw runtime_error(std::string("expected/") + realFilename + ".txt not found");

    vector<double> shortestLength;
    double i;
    while (file >> i) {
        shortestLength.push_back(i);
    }

    vector<double> lengths(shortestLength.size(), INFINITY);
    lengths[startIndex] = 0;
    for (const Funnel *const f : list) {
        if (lengths[f->p->index] > f->sp) lengths[f->p->index] = f->sp;
    }

    const double epsilon = 1e-9;
    const bool result = equal(shortestLength.begin(), shortestLength.end(), lengths.begin(),
                              [epsilon](const double a, const double b) { return fabs(a - b) < epsilon; });
    if (!result) {
        filesystem::create_directory("output");
        ofstream file(string("output/") + realFilename + " (" + to_string(n) + ").txt");
        file << fixed << setprecision(numeric_limits<double>::max_digits10);
        for (const double i : lengths) file << i << '\n';
    }

    return result;
}

void run(const char *filename, const size_t startIndex = 0) {
    cout << "File \"" << filename << "\": ";
    const pair<TriangleMesh, Point> meshAndStartPoint = getMeshAndStartPoint(filename, startIndex);

    const auto start = chrono::high_resolution_clock::now();
    const vector<Funnel*> list = FunnelTree(meshAndStartPoint.first, meshAndStartPoint.second);
    const auto end = chrono::high_resolution_clock::now();

    const auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Funnel tree with root " << startIndex << " initialized with " << list.size() << " nodes in " << duration.count() << " milliseconds.\n";
    if (!compareLength(filename, list, startIndex)) cout << "---------- NOT PASSED ----------\n";

    deleteFunnelTree(list);
}

void time(const char *filename, const size_t startIndex = 0, const size_t n = 100) {
    cout << "File \"" << filename << "\" ran " << n << " times. Avg: ";
    const pair<TriangleMesh, Point> meshAndStartPoint = getMeshAndStartPoint(filename, startIndex);

    vector<chrono::nanoseconds> durations;
    vector<bool> passed;
    for (size_t i = 0; i < n; i++) {
        const auto start = chrono::high_resolution_clock::now();
        const vector<Funnel*> list = FunnelTree(meshAndStartPoint.first, meshAndStartPoint.second);
        const auto end = chrono::high_resolution_clock::now();

        durations.emplace_back(end - start);
        passed.push_back(compareLength(filename, list, startIndex, i));
        deleteFunnelTree(list);
    }

    const chrono::nanoseconds avg = reduce(durations.begin(), durations.end()) / n;
    chrono::nanoseconds error(0);
    
    // rand error
    for (const chrono::nanoseconds i : durations) error += i < avg ? avg - i : i - avg;
    error /= n;

    error++;    // sys error
    cout << chrono::duration_cast<chrono::microseconds>(avg).count() << " +/- "
         << chrono::duration_cast<chrono::microseconds>(error).count() << " microseconds ("
         << count_if(passed.begin(), passed.end(), [](const bool i) { return i; }) << " passed)\n";
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