#pragma once
#include "grid.hpp"
void initialize_MHD_disk(FlowField& flow, int seed = 12345);
void add_divergence_error(FlowField& flow, double amplitude = 0.1); 
void initialize_orszag_tang(FlowField& flow);  // Add this line

