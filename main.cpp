#include "solver.hpp"
#include "physics.hpp"
#include "io.hpp"

#include <filesystem>
#include <chrono>
#include <iostream>
#include <iomanip>

static std::string prepare_output_dir(){
    namespace fs = std::filesystem;
    fs::path base("Result");
    if(fs::exists(base) && !fs::is_empty(base)){
        auto ts = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        fs::rename(base, "Result_"+std::to_string(ts));
    }
    fs::create_directory(base);
    return "Result";
}

int main(){
    const int nx=64, ny=64;
    const double Lx=1.0,Ly=1.0, dx=Lx/(nx-1), dy=Ly/(ny-1);
    const double nu=0.01;
    const int max_steps=2000;
    const int output_every=20;
    const double t_end = 20.0;

    std::string out_dir = prepare_output_dir();

    FlowField flow(nx,ny,dx,dy);
    initialize_orszag_tang(flow);
    //initialize_MHD_disk(flows[0]); // deterministic seed default
    //add_divergence_error(flows[0], 0.1);


    auto t0=std::chrono::high_resolution_clock::now();
    double t = 0.0;
    for(int step=0; step<=max_steps && t < t_end; ++step){
        // Use dynamic CFL-based timestep from the current flow state
        double dt = compute_cfl_timestep(flow);
        if(t + dt > t_end) dt = t_end - t;

        solve_MHD(flow, dt, nu);
        t += dt;

        if(step%output_every==0){
            auto [max_divB, L1_divB] = compute_divergence_errors(flow);
            std::cout << "step "<< std::setw(4) << step << " dt="<<dt
                      << " max_divB=" << max_divB
                      << " L1_divB=" << L1_divB << "\n";
            save_flow_MHD(flow,out_dir,step);
        }
    }
    auto t1=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;
    std::cout<<"Total time "<<elapsed.count()<<" s\n";
    return 0;
}
