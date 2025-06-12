#pragma once
#include <vector>
#include <cmath>
#include <string>

// Maximum number of AMR patches allowed in a single run
constexpr int MAX_REFINEMENTS = 5;
extern int refinement_count;

/**
 * Lightweight 2‑D uniformly‑spaced scalar field.
 */
class Grid {
public:
    int nx, ny;
    double dx, dy;
    double x0, y0;
    std::vector<std::vector<double>> data;

    Grid(int nx, int ny, double dx, double dy, double x0=0.0, double y0=0.0);
    void fill(double v);
};

class AMRGrid {
public:
    std::vector<Grid> levels;
    int max_level;
    AMRGrid(int base_nx, int base_ny, double Lx, double Ly, int max_level=1);
    void refine(int level, int start_x, int start_y, int fine_nx, int fine_ny);
    bool needs_refinement(const Grid& g, int i,int j,double threshold);
};

struct FlowField {
    Grid rho,u,v,p,e;
    Grid bx,by,psi;
    FlowField(int nx,int ny,double dx,double dy,double x0=0.0,double y0=0.0);
    FlowField(const Grid& g);
};
