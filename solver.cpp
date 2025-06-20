// new version
#include "solver.hpp"
#include <omp.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

static constexpr double ETA = 0.001;    // Magnetic diffusivity
static constexpr double CH = 0.8;      // GLM wave speed
static constexpr double CR = 0.01;     // GLM damping coefficient (improved value)
static constexpr double gamma_gas = 5.0/3.0;

// Helper function: compute Laplacian
static inline double laplacian(const Grid& g, int i, int j) {
    return (g.data[i+1][j] - 2*g.data[i][j] + g.data[i-1][j])/(g.dx*g.dx)
         + (g.data[i][j+1] - 2*g.data[i][j] + g.data[i][j-1])/(g.dy*g.dy);
}

// Minmod slope limiter
static inline double minmod(double a, double b){
    if(a*b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

// Compute fast magnetosonic speed (for CFL condition)
static double compute_fast_speed(double rho, double p, double Bx, double By) {
    double cs2 = gamma_gas * p / rho;  // Sound speed squared
    double ca2 = (Bx*Bx + By*By) / rho; // Alfven speed squared
    return sqrt(cs2 + ca2);
}

// HLL Riemann solver structure
struct HLLFlux {
    double F_rho, F_momx, F_momy, F_E, F_Bx, F_By, F_psi;
};

// HLL flux computation in X direction
HLLFlux compute_hll_flux_x(double rhoL, double uL, double vL, double pL, double BxL, double ByL, double psiL,
                           double rhoR, double uR, double vR, double pR, double BxR, double ByR, double psiR) {
    // Compute total pressure and energy
    double B2L = BxL*BxL + ByL*ByL;
    double B2R = BxR*BxR + ByR*ByR;
    double ptL = pL + 0.5*B2L;  // Total pressure
    double ptR = pR + 0.5*B2R;
    double EL = pL/(gamma_gas-1) + 0.5*rhoL*(uL*uL + vL*vL) + 0.5*B2L;
    double ER = pR/(gamma_gas-1) + 0.5*rhoR*(uR*uR + vR*vR) + 0.5*B2R;
    
    // Compute wave speeds
    double cfL = compute_fast_speed(rhoL, pL, BxL, ByL);
    double cfR = compute_fast_speed(rhoR, pR, BxR, ByR);
    double SL = std::min(uL - cfL, uR - cfR);
    double SR = std::max(uL + cfL, uR + cfR);
    
    HLLFlux flux;
    
    if (SL > 0) {
        // Left state flux
        flux.F_rho = rhoL * uL;
        flux.F_momx = rhoL * uL * uL + ptL - BxL * BxL;
        flux.F_momy = rhoL * uL * vL - BxL * ByL;
        flux.F_E = (EL + ptL) * uL - BxL * (uL*BxL + vL*ByL);
        flux.F_Bx = psiL;  // GLM
        flux.F_By = uL * ByL - vL * BxL;
        flux.F_psi = CH * CH * BxL;
    }
    else if (SR < 0) {
        // Right state flux
        flux.F_rho = rhoR * uR;
        flux.F_momx = rhoR * uR * uR + ptR - BxR * BxR;
        flux.F_momy = rhoR * uR * vR - BxR * ByR;
        flux.F_E = (ER + ptR) * uR - BxR * (uR*BxR + vR*ByR);
        flux.F_Bx = psiR;  // GLM
        flux.F_By = uR * ByR - vR * BxR;
        flux.F_psi = CH * CH * BxR;
    }
    else {
        // HLL average
        double FL_rho = rhoL * uL;
        double FR_rho = rhoR * uR;
        double FL_momx = rhoL * uL * uL + ptL - BxL * BxL;
        double FR_momx = rhoR * uR * uR + ptR - BxR * BxR;
        double FL_momy = rhoL * uL * vL - BxL * ByL;
        double FR_momy = rhoR * uR * vR - BxR * ByR;
        double FL_E = (EL + ptL) * uL - BxL * (uL*BxL + vL*ByL);
        double FR_E = (ER + ptR) * uR - BxR * (uR*BxR + vR*ByR);
        double FL_By = uL * ByL - vL * BxL;
        double FR_By = uR * ByR - vR * BxR;
        
        flux.F_rho = (SR * FL_rho - SL * FR_rho + SL * SR * (rhoR - rhoL)) / (SR - SL);
        flux.F_momx = (SR * FL_momx - SL * FR_momx + SL * SR * (rhoR*uR - rhoL*uL)) / (SR - SL);
        flux.F_momy = (SR * FL_momy - SL * FR_momy + SL * SR * (rhoR*vR - rhoL*vL)) / (SR - SL);
        flux.F_E = (SR * FL_E - SL * FR_E + SL * SR * (ER - EL)) / (SR - SL);
        flux.F_Bx = (SR * psiL - SL * psiR + SL * SR * (BxR - BxL)) / (SR - SL);
        flux.F_By = (SR * FL_By - SL * FR_By + SL * SR * (ByR - ByL)) / (SR - SL);
        flux.F_psi = CH * CH * (SR * BxL - SL * BxR + SL * SR * (psiR - psiL)) / (SR - SL);
    }
    
    return flux;
}

// HLL flux computation in Y direction (similar to X direction)
HLLFlux compute_hll_flux_y(double rhoL, double uL, double vL, double pL, double BxL, double ByL, double psiL,
                           double rhoR, double uR, double vR, double pR, double BxR, double ByR, double psiR) {
    double B2L = BxL*BxL + ByL*ByL;
    double B2R = BxR*BxR + ByR*ByR;
    double ptL = pL + 0.5*B2L;
    double ptR = pR + 0.5*B2R;
    double EL = pL/(gamma_gas-1) + 0.5*rhoL*(uL*uL + vL*vL) + 0.5*B2L;
    double ER = pR/(gamma_gas-1) + 0.5*rhoR*(uR*uR + vR*vR) + 0.5*B2R;
    
    double cfL = compute_fast_speed(rhoL, pL, BxL, ByL);
    double cfR = compute_fast_speed(rhoR, pR, BxR, ByR);
    double SL = std::min(vL - cfL, vR - cfR);
    double SR = std::max(vL + cfL, vR + cfR);
    
    HLLFlux flux;
    
    if (SL > 0) {
        flux.F_rho = rhoL * vL;
        flux.F_momx = rhoL * vL * uL - ByL * BxL;
        flux.F_momy = rhoL * vL * vL + ptL - ByL * ByL;
        flux.F_E = (EL + ptL) * vL - ByL * (uL*BxL + vL*ByL);
        flux.F_Bx = vL * BxL - uL * ByL;
        flux.F_By = psiL;  // GLM
        flux.F_psi = CH * CH * ByL;
    }
    else if (SR < 0) {
        flux.F_rho = rhoR * vR;
        flux.F_momx = rhoR * vR * uR - ByR * BxR;
        flux.F_momy = rhoR * vR * vR + ptR - ByR * ByR;
        flux.F_E = (ER + ptR) * vR - ByR * (uR*BxR + vR*ByR);
        flux.F_Bx = vR * BxR - uR * ByR;
        flux.F_By = psiR;  // GLM
        flux.F_psi = CH * CH * ByR;
    }
    else {
        // HLL average (similar to X direction)
        double FL_rho = rhoL * vL;
        double FR_rho = rhoR * vR;
        double FL_momx = rhoL * vL * uL - ByL * BxL;
        double FR_momx = rhoR * vR * uR - ByR * BxR;
        double FL_momy = rhoL * vL * vL + ptL - ByL * ByL;
        double FR_momy = rhoR * vR * vR + ptR - ByR * ByR;
        double FL_E = (EL + ptL) * vL - ByL * (uL*BxL + vL*ByL);
        double FR_E = (ER + ptR) * vR - ByR * (uR*BxR + vR*ByR);
        double FL_Bx = vL * BxL - uL * ByL;
        double FR_Bx = vR * BxR - uR * ByR;
        
        flux.F_rho = (SR * FL_rho - SL * FR_rho + SL * SR * (rhoR - rhoL)) / (SR - SL);
        flux.F_momx = (SR * FL_momx - SL * FR_momx + SL * SR * (rhoR*uR - rhoL*uL)) / (SR - SL);
        flux.F_momy = (SR * FL_momy - SL * FR_momy + SL * SR * (rhoR*vR - rhoL*vL)) / (SR - SL);
        flux.F_E = (SR * FL_E - SL * FR_E + SL * SR * (ER - EL)) / (SR - SL);
        flux.F_Bx = (SR * FL_Bx - SL * FR_Bx + SL * SR * (BxR - BxL)) / (SR - SL);
        flux.F_By = (SR * psiL - SL * psiR + SL * SR * (ByR - ByL)) / (SR - SL);
        flux.F_psi = CH * CH * (SR * ByL - SL * ByR + SL * SR * (psiR - psiL)) / (SR - SL);
    }
    
    return flux;
}

// Compute dynamic CFL timestep
double compute_cfl_timestep(const FlowField& flow, double cfl_number) {
    double dt_min = 1e10;
    const Grid& grid = flow.rho;
    
    #pragma omp parallel for collapse(2) reduction(min:dt_min)
    for (int i = 1; i < grid.nx-1; ++i) {
        for (int j = 1; j < grid.ny-1; ++j) {
            double rho = flow.rho.data[i][j];
            double u = flow.u.data[i][j];
            double v = flow.v.data[i][j];
            double p = flow.p.data[i][j];
            double Bx = flow.bx.data[i][j];
            double By = flow.by.data[i][j];
            
            double cf = compute_fast_speed(rho, p, Bx, By);
            
            double dt_x = grid.dx / (std::abs(u) + cf);
            double dt_y = grid.dy / (std::abs(v) + cf);
            
            dt_min = std::min(dt_min, std::min(dt_x, dt_y));
        }
    }
    
    double dt_glm = std::min(grid.dx, grid.dy) / CH;
    if(dt_min > 1.0) // prevent unrealistically large dt due to NaNs
        dt_min = std::min(grid.dx, grid.dy) / CH;
    return cfl_number * std::min(dt_min, dt_glm);
}

// Compute divergence errors for monitoring
std::pair<double, double> compute_divergence_errors(const FlowField& flow) {
    const Grid& grid = flow.bx;
    double max_divB = 0.0;
    double L1_divB = 0.0;
    int count = 0;
    
    #pragma omp parallel for collapse(2) reduction(max:max_divB) reduction(+:L1_divB,count)
    for (int i = 1; i < grid.nx-1; ++i) {
        for (int j = 1; j < grid.ny-1; ++j) {
            double divB = (flow.bx.data[i+1][j] - flow.bx.data[i-1][j]) / (2*grid.dx)
                        + (flow.by.data[i][j+1] - flow.by.data[i][j-1]) / (2*grid.dy);
            
            double abs_divB = std::abs(divB);
            max_divB = std::max(max_divB, abs_divB);
            L1_divB += abs_divB;
            count++;
        }
    }
    
    L1_divB /= count;
    return {max_divB, L1_divB};
}

// Main improved MHD solver function

static void update_level(FlowField& flow,double dt,double nu){
    Grid& grid = flow.rho;
    
    // Use dynamic CFL timestep
    double dt_cfl = compute_cfl_timestep(flow);
    dt = std::min(dt, dt_cfl);
    
    // Temporary arrays
    auto rho_new = flow.rho.data;
    auto momx_new = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto momy_new = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto e_new = flow.e.data;
    auto bx_new = flow.bx.data;
    auto by_new = flow.by.data;
    auto psi_new = flow.psi.data;

    // Pre-compute limited slopes for MUSCL reconstruction
    auto srho_x = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto su_x   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sv_x   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sp_x   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sbx_x  = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sby_x  = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto spsi_x = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));

    auto srho_y = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto su_y   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sv_y   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sp_y   = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sbx_y  = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto sby_y  = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));
    auto spsi_y = std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.ny));

    auto zero_slopes = [](auto& arr){
        for(auto& r : arr) std::fill(r.begin(), r.end(), 0.0);
    };
    zero_slopes(srho_x);  zero_slopes(su_x);  zero_slopes(sv_x);
    zero_slopes(sp_x);    zero_slopes(sbx_x); zero_slopes(sby_x); zero_slopes(spsi_x);
    zero_slopes(srho_y);  zero_slopes(su_y);  zero_slopes(sv_y);
    zero_slopes(sp_y);    zero_slopes(sbx_y); zero_slopes(sby_y); zero_slopes(spsi_y);

    // Slopes in X direction
    #pragma omp parallel for collapse(2)
    for(int i=1;i<grid.nx-1;++i){
        for(int j=0;j<grid.ny;++j){
            srho_x[i][j] = minmod(flow.rho.data[i][j]-flow.rho.data[i-1][j],
                                  flow.rho.data[i+1][j]-flow.rho.data[i][j]);
            su_x[i][j]   = minmod(flow.u.data[i][j]-flow.u.data[i-1][j],
                                  flow.u.data[i+1][j]-flow.u.data[i][j]);
            sv_x[i][j]   = minmod(flow.v.data[i][j]-flow.v.data[i-1][j],
                                  flow.v.data[i+1][j]-flow.v.data[i][j]);
            sp_x[i][j]   = minmod(flow.p.data[i][j]-flow.p.data[i-1][j],
                                  flow.p.data[i+1][j]-flow.p.data[i][j]);
            sbx_x[i][j]  = minmod(flow.bx.data[i][j]-flow.bx.data[i-1][j],
                                  flow.bx.data[i+1][j]-flow.bx.data[i][j]);
            sby_x[i][j]  = minmod(flow.by.data[i][j]-flow.by.data[i-1][j],
                                  flow.by.data[i+1][j]-flow.by.data[i][j]);
            spsi_x[i][j] = minmod(flow.psi.data[i][j]-flow.psi.data[i-1][j],
                                  flow.psi.data[i+1][j]-flow.psi.data[i][j]);
        }
    }

    // Slopes in Y direction
    #pragma omp parallel for collapse(2)
    for(int i=0;i<grid.nx;++i){
        for(int j=1;j<grid.ny-1;++j){
            srho_y[i][j] = minmod(flow.rho.data[i][j]-flow.rho.data[i][j-1],
                                  flow.rho.data[i][j+1]-flow.rho.data[i][j]);
            su_y[i][j]   = minmod(flow.u.data[i][j]-flow.u.data[i][j-1],
                                  flow.u.data[i][j+1]-flow.u.data[i][j]);
            sv_y[i][j]   = minmod(flow.v.data[i][j]-flow.v.data[i][j-1],
                                  flow.v.data[i][j+1]-flow.v.data[i][j]);
            sp_y[i][j]   = minmod(flow.p.data[i][j]-flow.p.data[i][j-1],
                                  flow.p.data[i][j+1]-flow.p.data[i][j]);
            sbx_y[i][j]  = minmod(flow.bx.data[i][j]-flow.bx.data[i][j-1],
                                  flow.bx.data[i][j+1]-flow.bx.data[i][j]);
            sby_y[i][j]  = minmod(flow.by.data[i][j]-flow.by.data[i][j-1],
                                  flow.by.data[i][j+1]-flow.by.data[i][j]);
            spsi_y[i][j] = minmod(flow.psi.data[i][j]-flow.psi.data[i][j-1],
                                  flow.psi.data[i][j+1]-flow.psi.data[i][j]);
        }
    }
    
    // First compute momentum (for HLL solver)
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            momx_new[i][j] = flow.rho.data[i][j] * flow.u.data[i][j];
            momy_new[i][j] = flow.rho.data[i][j] * flow.v.data[i][j];
        }
    }
    
    // Update using HLL solver
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < grid.nx-1; ++i) {
        for (int j = 1; j < grid.ny-1; ++j) {
            // Get current state
            double rho = flow.rho.data[i][j];
            double u = flow.u.data[i][j];
            double v = flow.v.data[i][j];
            double p = flow.p.data[i][j];
            double Bx = flow.bx.data[i][j];
            double By = flow.by.data[i][j];
            double psi = flow.psi.data[i][j];
            
            // X direction fluxes with MUSCL reconstruction
            HLLFlux flux_xp = compute_hll_flux_x(
               // left state at i+1/2
                rho + 0.5*srho_x[i][j],
                u   + 0.5*su_x[i][j],
                v   + 0.5*sv_x[i][j],
                p   + 0.5*sp_x[i][j],
                Bx  + 0.5*sbx_x[i][j],
                By  + 0.5*sby_x[i][j],
                psi + 0.5*spsi_x[i][j],
                // right state at i+1/2
                flow.rho.data[i+1][j] - 0.5*srho_x[i+1][j],
                flow.u.data[i+1][j]   - 0.5*su_x[i+1][j],
                flow.v.data[i+1][j]   - 0.5*sv_x[i+1][j],
                flow.p.data[i+1][j]   - 0.5*sp_x[i+1][j],
                flow.bx.data[i+1][j]  - 0.5*sbx_x[i+1][j],
                flow.by.data[i+1][j]  - 0.5*sby_x[i+1][j],
                flow.psi.data[i+1][j] - 0.5*spsi_x[i+1][j]
            );
            
            HLLFlux flux_xm = compute_hll_flux_x(
                // left state at i-1/2
                flow.rho.data[i-1][j] + 0.5*srho_x[i-1][j],
                flow.u.data[i-1][j]   + 0.5*su_x[i-1][j],
                flow.v.data[i-1][j]   + 0.5*sv_x[i-1][j],
                flow.p.data[i-1][j]   + 0.5*sp_x[i-1][j],
                flow.bx.data[i-1][j]  + 0.5*sbx_x[i-1][j],
                flow.by.data[i-1][j]  + 0.5*sby_x[i-1][j],
                flow.psi.data[i-1][j] + 0.5*spsi_x[i-1][j],
                // right state at i-1/2
                rho - 0.5*srho_x[i][j],
                u   - 0.5*su_x[i][j],
                v   - 0.5*sv_x[i][j],
                p   - 0.5*sp_x[i][j],
                Bx  - 0.5*sbx_x[i][j],
                By  - 0.5*sby_x[i][j],
                psi - 0.5*spsi_x[i][j]
            );
            
            // Y direction fluxes with MUSCL reconstruction
            HLLFlux flux_yp = compute_hll_flux_y(
                // bottom state at j+1/2
                rho + 0.5*srho_y[i][j],
                u   + 0.5*su_y[i][j],
                v   + 0.5*sv_y[i][j],
                p   + 0.5*sp_y[i][j],
                Bx  + 0.5*sbx_y[i][j],
                By  + 0.5*sby_y[i][j],
                psi + 0.5*spsi_y[i][j],
                // top state at j+1/2
                flow.rho.data[i][j+1] - 0.5*srho_y[i][j+1],
                flow.u.data[i][j+1]   - 0.5*su_y[i][j+1],
                flow.v.data[i][j+1]   - 0.5*sv_y[i][j+1],
                flow.p.data[i][j+1]   - 0.5*sp_y[i][j+1],
                flow.bx.data[i][j+1]  - 0.5*sbx_y[i][j+1],
                flow.by.data[i][j+1]  - 0.5*sby_y[i][j+1],
                flow.psi.data[i][j+1] - 0.5*spsi_y[i][j+1]
            );
            
            HLLFlux flux_ym = compute_hll_flux_y(
                // bottom state at j-1/2
                flow.rho.data[i][j-1] + 0.5*srho_y[i][j-1],
                flow.u.data[i][j-1]   + 0.5*su_y[i][j-1],
                flow.v.data[i][j-1]   + 0.5*sv_y[i][j-1],
                flow.p.data[i][j-1]   + 0.5*sp_y[i][j-1],
                flow.bx.data[i][j-1]  + 0.5*sbx_y[i][j-1],
                flow.by.data[i][j-1]  + 0.5*sby_y[i][j-1],
                flow.psi.data[i][j-1] + 0.5*spsi_y[i][j-1],
                // top state at j-1/2
                rho - 0.5*srho_y[i][j],
                u   - 0.5*su_y[i][j],
                v   - 0.5*sv_y[i][j],
                p   - 0.5*sp_y[i][j],
                Bx  - 0.5*sbx_y[i][j],
                By  - 0.5*sby_y[i][j],
                psi - 0.5*spsi_y[i][j]
            );
            
            // Update conserved variables
            rho_new[i][j] = rho - dt/grid.dx * (flux_xp.F_rho - flux_xm.F_rho)
                                - dt/grid.dy * (flux_yp.F_rho - flux_ym.F_rho);
            
            momx_new[i][j] = momx_new[i][j] - dt/grid.dx * (flux_xp.F_momx - flux_xm.F_momx)
                                             - dt/grid.dy * (flux_yp.F_momx - flux_ym.F_momx);
            
            momy_new[i][j] = momy_new[i][j] - dt/grid.dx * (flux_xp.F_momy - flux_xm.F_momy)
                                             - dt/grid.dy * (flux_yp.F_momy - flux_ym.F_momy);
            
            e_new[i][j] = flow.e.data[i][j] - dt/grid.dx * (flux_xp.F_E - flux_xm.F_E)
                                            - dt/grid.dy * (flux_yp.F_E - flux_ym.F_E);
            double ke_temp = 0.5 * rho_new[i][j] * (u*u + v*v);
            double me_temp = 0.5 * (bx_new[i][j]*bx_new[i][j] + by_new[i][j]*by_new[i][j]);
            if (e_new[i][j] < ke_temp + me_temp + 1e-10) {
                std::cerr << "Warning: Insufficient total energy at ("<<i<<","<<j<<"), adjusting\n";
                e_new[i][j] = ke_temp + me_temp + 1e-10;
            }
            
            bx_new[i][j] = Bx - dt/grid.dx * (flux_xp.F_Bx - flux_xm.F_Bx)
                              - dt/grid.dy * (flux_yp.F_Bx - flux_ym.F_Bx);
            
            by_new[i][j] = By - dt/grid.dx * (flux_xp.F_By - flux_xm.F_By)
                              - dt/grid.dy * (flux_yp.F_By - flux_ym.F_By);
            
            psi_new[i][j] = psi - dt/grid.dx * (flux_xp.F_psi - flux_xm.F_psi)
                                - dt/grid.dy * (flux_yp.F_psi - flux_ym.F_psi);
            
            // Add viscous terms
            if (nu > 0) {
                momx_new[i][j] += dt * nu * rho * laplacian(flow.u, i, j);
                momy_new[i][j] += dt * nu * rho * laplacian(flow.v, i, j);
            }
            
            // Add magnetic diffusion
            if (ETA > 0) {
                bx_new[i][j] += dt * ETA * laplacian(flow.bx, i, j);
                by_new[i][j] += dt * ETA * laplacian(flow.by, i, j);
            }
            
            // GLM flux part handled above; divergence cleaning will be applied later
            
            // Ensure physical values
            rho_new[i][j] = std::max(rho_new[i][j], 1e-10);
            e_new[i][j]   = std::max(e_new[i][j], 1e-10);
        }
    }
    
    // Update primitive variables
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < grid.nx-1; ++i) {
        for (int j = 1; j < grid.ny-1; ++j) {
            flow.rho.data[i][j] = rho_new[i][j];
            flow.u.data[i][j] = momx_new[i][j] / rho_new[i][j];
            flow.v.data[i][j] = momy_new[i][j] / rho_new[i][j];
            flow.bx.data[i][j] = bx_new[i][j];
            flow.by.data[i][j] = by_new[i][j];
            flow.e.data[i][j]  = e_new[i][j];
            
            // Update pressure
            double ke = 0.5 * rho_new[i][j] * (flow.u.data[i][j]*flow.u.data[i][j] +
                                                flow.v.data[i][j]*flow.v.data[i][j]);
            double me = 0.5 * (bx_new[i][j]*bx_new[i][j] + by_new[i][j]*by_new[i][j]);
            double ie = e_new[i][j] - ke - me;
            if (ie < 0)
                std::cerr << "Warning: Negative internal energy at ("<<i<<","<<j<<")\n";
            flow.p.data[i][j] = (gamma_gas - 1.0) * std::max(ie, 1e-10);
        }
    }
    
    // Boundary conditions (periodic)
    #pragma omp parallel for
    for (int j = 0; j < grid.ny; ++j) {
        // X direction periodic BC using modulo indices
        int left_src  = (grid.nx + 0 - 2) % grid.nx;  // nx-2
        int right_src = (grid.nx - 1 + 2) % grid.nx;  // 1
        flow.rho.data[0][j] = flow.rho.data[left_src][j];
        flow.rho.data[grid.nx-1][j] = flow.rho.data[right_src][j];
        flow.u.data[0][j] = flow.u.data[left_src][j];
        flow.u.data[grid.nx-1][j] = flow.u.data[right_src][j];
        flow.v.data[0][j] = flow.v.data[left_src][j];
        flow.v.data[grid.nx-1][j] = flow.v.data[right_src][j];
        flow.p.data[0][j] = flow.p.data[left_src][j];
        flow.p.data[grid.nx-1][j] = flow.p.data[right_src][j];
        flow.e.data[0][j] = flow.e.data[left_src][j];
        flow.e.data[grid.nx-1][j] = flow.e.data[right_src][j];
        flow.bx.data[0][j] = flow.bx.data[left_src][j];
        flow.bx.data[grid.nx-1][j] = flow.bx.data[right_src][j];
        flow.by.data[0][j] = flow.by.data[left_src][j];
        flow.by.data[grid.nx-1][j] = flow.by.data[right_src][j];
        flow.psi.data[0][j] = flow.psi.data[left_src][j];
        flow.psi.data[grid.nx-1][j] = flow.psi.data[right_src][j];
    }
    
    #pragma omp parallel for
    for (int i = 0; i < grid.nx; ++i) {
        // Y direction periodic BC using modulo indices
        int bot_src = (grid.ny + 0 - 2) % grid.ny;  // ny-2
        int top_src = (grid.ny - 1 + 2) % grid.ny;  // 1
        flow.rho.data[i][0] = flow.rho.data[i][bot_src];
        flow.rho.data[i][grid.ny-1] = flow.rho.data[i][top_src];
        flow.u.data[i][0] = flow.u.data[i][bot_src];
        flow.u.data[i][grid.ny-1] = flow.u.data[i][top_src];
        flow.v.data[i][0] = flow.v.data[i][bot_src];
        flow.v.data[i][grid.ny-1] = flow.v.data[i][top_src];
        flow.p.data[i][0] = flow.p.data[i][bot_src];
        flow.p.data[i][grid.ny-1] = flow.p.data[i][top_src];
        flow.e.data[i][0] = flow.e.data[i][bot_src];
        flow.e.data[i][grid.ny-1] = flow.e.data[i][top_src];
        flow.bx.data[i][0] = flow.bx.data[i][bot_src];
        flow.bx.data[i][grid.ny-1] = flow.bx.data[i][top_src];
        flow.by.data[i][0] = flow.by.data[i][bot_src];
        flow.by.data[i][grid.ny-1] = flow.by.data[i][top_src];
        flow.psi.data[i][0] = flow.psi.data[i][bot_src];
        flow.psi.data[i][grid.ny-1] = flow.psi.data[i][top_src];
    }

    // GLM divergence cleaning including boundaries
    #pragma omp parallel for collapse(2)
    for(int i=0;i<grid.nx;++i){
        for(int j=0;j<grid.ny;++j){
            int ip=(i+1)%grid.nx, im=(i-1+grid.nx)%grid.nx;
            int jp=(j+1)%grid.ny, jm=(j-1+grid.ny)%grid.ny;
            double divB_new = (bx_new[ip][j] - bx_new[im][j])/(2*grid.dx)
                            + (by_new[i][jp] - by_new[i][jm])/(2*grid.dy);
            flow.psi.data[i][j] = psi_new[i][j] - dt*CH*CH*divB_new
                                   - dt*CR*psi_new[i][j];
        }
    }
}

void solve_MHD(FlowField& flow, double dt, double nu){
    update_level(flow, dt, nu);
}
