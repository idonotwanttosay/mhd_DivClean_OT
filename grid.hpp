#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>


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


struct FlowField {
    Grid rho,u,v,p,e;
    Grid bx,by,psi;
    FlowField(int nx,int ny,double dx,double dy,double x0=0.0,double y0=0.0);
    FlowField(const Grid& g);
};
