#pragma once
#include "grid.hpp"
#include <vector>
void solve_MHD(AMRGrid& amr, std::vector<FlowField>& flows,double dt,double nu,int max_iter,double tol);
// Estimate stable timestep based on CFL condition
double compute_cfl_timestep(const FlowField& flow, double cfl_number = 0.4);
std::pair<double, double> compute_divergence_errors(const FlowField& flow);  // Add this line
