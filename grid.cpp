#include "grid.hpp"
#include <algorithm>
#include <omp.h>
#include <iostream>

Grid::Grid(int nx_,int ny_,double dx_,double dy_,double x0_,double y0_)
    : nx(nx_),ny(ny_),dx(dx_),dy(dy_),x0(x0_),y0(y0_),
      data(nx_, std::vector<double>(ny_,0.0))
{
    if(nx_ < 3 || ny_ < 3)
        throw std::invalid_argument("Grid size must be at least 3x3");
}

void Grid::fill(double v){
#pragma omp parallel for collapse(2)
    for(int i=0;i<nx;++i)
        for(int j=0;j<ny;++j)
            data[i][j]=v;
}


FlowField::FlowField(int nx,int ny,double dx,double dy,double x0,double y0)
    : rho(nx,ny,dx,dy,x0,y0), u(nx,ny,dx,dy,x0,y0), v(nx,ny,dx,dy,x0,y0),
      p(nx,ny,dx,dy,x0,y0), e(nx,ny,dx,dy,x0,y0),
      bx(nx,ny,dx,dy,x0,y0), by(nx,ny,dx,dy,x0,y0), psi(nx,ny,dx,dy,x0,y0)
{
    if(nx < 3 || ny < 3)
        throw std::invalid_argument("FlowField grid must be at least 3x3");
}

FlowField::FlowField(const Grid& g)
    : FlowField(g.nx,g.ny,g.dx,g.dy,g.x0,g.y0) {}
