#include "grid.hpp"
#include <algorithm>
#include <omp.h>
#include <iostream>

int refinement_count = 0;

Grid::Grid(int nx_,int ny_,double dx_,double dy_,double x0_,double y0_)
    : nx(nx_),ny(ny_),dx(dx_),dy(dy_),x0(x0_),y0(y0_),
      data(nx_, std::vector<double>(ny_,0.0)) {}

void Grid::fill(double v){
#pragma omp parallel for collapse(2)
    for(int i=0;i<nx;++i)
        for(int j=0;j<ny;++j)
            data[i][j]=v;
}

AMRGrid::AMRGrid(int base_nx,int base_ny,double Lx,double Ly,int max_level_)
    : max_level(max_level_)
{
    levels.emplace_back(base_nx,base_ny,Lx/(base_nx-1),Ly/(base_ny-1),0.0,0.0);
}

void AMRGrid::refine(int level, int start_x, int start_y, int coarse_nx_patch, int coarse_ny_patch) {
    if(level >= max_level || refinement_count >= MAX_REFINEMENTS) return;
    
    // Ensure the patch stays within bounds
    start_x = std::max(0, std::min(start_x, levels[level].nx - coarse_nx_patch));
    start_y = std::max(0, std::min(start_y, levels[level].ny - coarse_ny_patch));
    
    // Fine grid has 2x resolution
    int fine_nx = coarse_nx_patch * 2 + 1;  // +1 for proper overlap
    int fine_ny = coarse_ny_patch * 2 + 1;
    
    double fine_dx = levels[level].dx / 2.0;
    double fine_dy = levels[level].dy / 2.0;
    
    // Align with coarse cell edges
    double x0 = levels[level].x0 + start_x * levels[level].dx;
    double y0 = levels[level].y0 + start_y * levels[level].dy;
    
    std::cout << "[AMR] refine level " << level << " at (" << start_x << "," << start_y 
              << ") covering " << coarse_nx_patch << "x" << coarse_ny_patch << " coarse cells\n";
    
    levels.emplace_back(fine_nx, fine_ny, fine_dx, fine_dy, x0, y0);
    ++refinement_count;
    
    // Proper interpolation
    Grid& fine = levels.back();
    Grid& coarse = levels[level];
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < fine_nx; ++i) {
        for(int j = 0; j < fine_ny; ++j) {
            // Map fine grid point to coarse grid
            double x_fine = x0 + i * fine_dx;
            double y_fine = y0 + j * fine_dy;
            
            // Find coarse cell indices
            double ci_real = (x_fine - coarse.x0) / coarse.dx;
            double cj_real = (y_fine - coarse.y0) / coarse.dy;
            
            int ci = (int)ci_real;
            int cj = (int)cj_real;
            
            // Bilinear interpolation
            if(ci >= 0 && ci < coarse.nx-1 && cj >= 0 && cj < coarse.ny-1) {
                double fx = ci_real - ci;
                double fy = cj_real - cj;
                
                fine.data[i][j] = (1-fx)*(1-fy)*coarse.data[ci][cj] +
                                  fx*(1-fy)*coarse.data[ci+1][cj] +
                                  (1-fx)*fy*coarse.data[ci][cj+1] +
                                  fx*fy*coarse.data[ci+1][cj+1];
            } else {
                // Boundary handling
                ci = std::max(0, std::min(ci, coarse.nx-1));
                cj = std::max(0, std::min(cj, coarse.ny-1));
                fine.data[i][j] = coarse.data[ci][cj];
            }
        }
    }
}

bool AMRGrid::needs_refinement(const Grid& g,int i,int j,double thres){
    if(i<=0||i>=g.nx-1||j<=0||j>=g.ny-1) return false;
    double gx=std::fabs(g.data[i+1][j]-g.data[i-1][j])/(2*g.dx);
    double gy=std::fabs(g.data[i][j+1]-g.data[i][j-1])/(2*g.dy);
    return (gx+gy)>thres;
}

FlowField::FlowField(int nx,int ny,double dx,double dy,double x0,double y0)
    : rho(nx,ny,dx,dy,x0,y0), u(nx,ny,dx,dy,x0,y0), v(nx,ny,dx,dy,x0,y0),
      p(nx,ny,dx,dy,x0,y0), e(nx,ny,dx,dy,x0,y0),
      bx(nx,ny,dx,dy,x0,y0), by(nx,ny,dx,dy,x0,y0), psi(nx,ny,dx,dy,x0,y0) {}

FlowField::FlowField(const Grid& g)
    : FlowField(g.nx,g.ny,g.dx,g.dy,g.x0,g.y0) {}
