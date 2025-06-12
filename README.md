# mhd_DivClean_OT

This repository contains a simple magnetohydrodynamics solver with OpenMP acceleration.

To build the executable, run the provided `compile.sh` script or use the following command:

```bash
g++ main.cpp grid.cpp physics.cpp solver.cpp io.cpp -std=c++17 -O2 -fopenmp -o mhd_solver
```

To run the solver and generate analysis plots, execute:

```bash
bash run.sh
```

The script builds `mhd_solver` if needed and runs the simulation. Output CSV
files are written to a new `Result` directory. If a previous `Result` folder
exists it is renamed with a timestamp. After validating the output,
`analysis_summary.py` generates summary plots and `plot_flow.py` creates an
animation of the flow field. All console output is stored in `solver.log`.
