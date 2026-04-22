#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    amrex::Print() << "Dim= " << AMREX_SPACEDIM << '\n';
    
    int max_step = -1;
    Real strt_time = 0.0;
    Real stop_time = 0.1;
    
    {   
        ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("strt_time", strt_time);
        pp.query("stop_time", stop_time);
    }
    
    amrex::Print() << "max_step: " << max_step << '\n';
    amrex::Print() << "strt_time: " << strt_time << '\n';
    amrex::Print() << "stop_time: " << stop_time << '\n';
    
    amrex::Finalize();
    return 0;
}