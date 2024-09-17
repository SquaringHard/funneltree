// only compile this when results in expected folder are incorrect or missing and you need to re-run getExpectedLength.exe

// #define I_AM_READY_TO_COMPILE_THIS_FILE

#ifdef I_AM_READY_TO_COMPILE_THIS_FILE
#include "ft.h"
#include <utility>      // pair
#include <string>
#include <fstream>      // ifstream, ofstream
#include <cmath>        // INFINITY
#include <filesystem>
#include <iostream>     // cin, cout
#include <sstream>      // istringstream
#include <iterator>     // istream_iterator
#include <iomanip>      // fixed, setprecision
#include <limits>       // numeric_limits

int main() {
    cout << "Enter the filenames (in folder \"input\"), separated by spaces, or press Enter to run all files: ";
    string filenames;
    getline(cin, filenames);

    istringstream iss(filenames);
    vector<string> filenamesVec(istream_iterator<string>{iss}, istream_iterator<string>{});
    if (filenamesVec.empty()) for (const auto &entry : filesystem::directory_iterator("input")) filenamesVec.push_back(entry.path().filename().string());

    const size_t startIndex = 0;
    for (const string& filename : filenamesVec) {
        const string realFilename = string(filename) + "_s=" + to_string(startIndex) + ".txt"; 
        cout << "\nMaking " << realFilename << " ...";

        ifstream input(string("input/") + filename);
        size_t v, f, e;
        input >> v >> f >> e;

        vector<Point> points;
        for (size_t i = 0; i < v; i++) {
            double x, y, z;
            input >> x >> y >> z;
            points.push_back({x, y, z});
        }

        vector<Triangle> trianglesPointsIndexes;
        for (size_t i = 0; i < f; i++) {
            indexType a, b, c;
            short three;
            input >> three >> a >> b >> c;
            trianglesPointsIndexes.push_back({a, b, c});
        }

        const vector<Funnel> list = FunnelTree(TriangleMesh(points, trianglesPointsIndexes), startIndex);

        vector<double> shortestLength(v, INFINITY);
        shortestLength[startIndex] = 0;
        for (const Funnel &f : list) if (shortestLength[f.p] > f.sp) shortestLength[f.p] = f.sp;

        filesystem::create_directory("expected");
        ofstream output(string("expected/") + realFilename);
        output << fixed << setprecision(numeric_limits<double>::max_digits10);
        for (const double i : shortestLength) output << i << '\n';
    }
}
#endif