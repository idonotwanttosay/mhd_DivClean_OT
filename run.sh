#!/bin/bash
# Compile and run the solver, then execute Python analysis scripts.
set -e

# Build solver if needed
if [ ! -x mhd_solver ]; then
    bash compile.sh
fi

# Run solver and capture output
./mhd_solver | tee solver.log

# Verify that CSV output files exist
if ! ls Result/out_*.csv >/dev/null 2>&1; then
    echo "No output generated in Result/" >&2
    exit 1
fi

# Run analysis scripts
python3 analysis_summary.py | tee -a solver.log
python3 plot_flow.py | tee -a solver.log
echo "Results stored in Result/ (previous runs are timestamped)."
