#include "io.hpp"
#include <fstream>
#include <filesystem>

static void dump_scalar(const Grid& g,const std::string& fname){
    std::ofstream out(fname);
    for(int i=0;i<g.nx;++i)
        for(int j=0;j<g.ny;++j){
            double x=g.x0+i*g.dx;
            double y=g.y0+j*g.dy;
            out<<x<<','<<y<<','<<g.data[i][j]<<'\n';
        }
}

void save_flow_MHD(const FlowField& flow,const std::string& dir,int step){
    std::filesystem::create_directory(dir);
    const std::string prefix = dir + "/out_";
    dump_scalar(flow.rho, prefix+"rho_"+std::to_string(step)+".csv");
    dump_scalar(flow.u,   prefix+"u_"+std::to_string(step)+".csv");
    dump_scalar(flow.v,   prefix+"v_"+std::to_string(step)+".csv");
    dump_scalar(flow.e,   prefix+"e_"+std::to_string(step)+".csv");
    dump_scalar(flow.bx,  prefix+"bx_"+std::to_string(step)+".csv");
    dump_scalar(flow.by,  prefix+"by_"+std::to_string(step)+".csv");
}
