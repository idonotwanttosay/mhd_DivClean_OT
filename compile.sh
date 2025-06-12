#!/bin/bash
# Compile the magnetohydrodynamics solver

g++ main.cpp grid.cpp physics.cpp solver.cpp io.cpp -std=c++17 -O2 -fopenmp -o mhd_solver
