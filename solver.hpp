#pragma once
#include "grid.hpp"
void solve_MHD(FlowField& flow, double dt, double nu);
// Estimate stable timestep based on CFL condition
double compute_cfl_timestep(const FlowField& flow, double cfl_number = 0.2);
std::pair<double, double> compute_divergence_errors(const FlowField& flow);  // Add this line
