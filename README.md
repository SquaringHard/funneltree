# FunnelTree algorithm in C++ and parallelized
The funnel tree algorithm is an algorithm for finding shortest paths on polyhedral surfaces. The article about the algorithm is at https://doi.org/10.1080/02331934.2023.2241496

# How to use
- compile files: `g++ -o demo main.cpp ft.cpp -fopenmp` (remove `fopenmp` for serial version)
- run file: `.\demo`, `.\demo spiral1.geom`, `.\demo cube1.geom spiral2.geom sphere4.geom`
- add datasets in input (refer to the already included .geom files for formats)
- get lengths of shortest path: run getExpectedLength.exe
- compare length: use `#define LENGTH_COMPARE` in main.cpp, run getExpectedLength.exe first
- get each thread's runtime: use `#define THREAD_TIMING` in ft.h
- run dataset in 100 times: use `time()` instead of `run()` in `main` function in main.cpp

# Contribution
[My friend](https://github.com/BanAnA9205)
