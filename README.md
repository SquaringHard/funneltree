# FunnelTree algorithm in C++ and parallelized
The funnel tree algorithm is an algorithm for finding shortest paths on polyhedral surfaces. The article about the algorithm is at https://doi.org/10.1080/02331934.2023.2241496

# How to use demo
- compile files: `g++ -o demo main.cpp ft.cpp -fopenmp` (remove `fopenmp` for serial version)
- run file: `.\demo`, `.\demo spiral1.geom`, `.\demo cube1.geom spiral2.geom sphere4.geom`
- add datasets in input folder (check the already included .geom files for formats)
- get lengths of shortest path: compile and run `getExpectedLength.cpp` instead of `main.cpp`
- compare length: use `#define LENGTH_COMPARE` in main.cpp, run getExpectedLength.exe first
- run dataset 100 times: replace `run` with `repeat` in `main` function in main.cpp

# Contribution
[My friend](https://github.com/BanAnA9205)
