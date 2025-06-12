# mhd_DivClean_OT

This repository contains a simple magnetohydrodynamics solver with OpenMP acceleration.

To build the executable, run the provided `compile.sh` script or use the following command:

```bash
g++ main.cpp grid.cpp physics.cpp solver.cpp io.cpp -std=c++17 -O2 -fopenmp -o mhd_solver
```
