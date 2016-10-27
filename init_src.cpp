#include <rsf.hh>
#include <cassert>
#include "grid.hpp"
#include "init.hpp"
#include <iostream>

void InitSrc(//output
             int& src_grid,
             int src_loc[],
             std::valarray<float>& src_trace,
             //input
             Grid grids[],
             int num_grids,
             real dt)
{
    iRSF src("src");
    
    //get source data
    int src_nt;
    real src_dt;
    src.get("n1",src_nt,0);
    src.get("d1",src_dt,0.); //unit [s]
    
    src_trace.resize(src_nt);
    src >> src_trace;
    
    //get source location
    real loc[3];
    src.get("o2",loc[0]);
    src.get("o3",loc[1]);
    src.get("o4",loc[2]);
    
    //find out where on grid we should put this source
    src_grid = 0;
    real tmp = loc[2];
    while(tmp > grids[src_grid].step[2]*grids[src_grid].size[2])
    {
        tmp -= grids[src_grid].step[2]*grids[src_grid].size[2];
        src_grid++;
    }
    assert(src_grid < num_grids);
    
    for(int i = 0; i < 2; i++)
        src_loc[i] = round((loc[i]-grids[src_grid].origin[i])/grids[src_grid].step[i]);
    src_loc[2] = round(tmp/grids[src_grid].step[2]);
    
    //print source location
    std::cout << "== Source is placed at... " << std::endl;
    std::cout << "(grid#" << src_grid << ", ix=" << src_loc[0] << ", iy=" << src_loc[1] << ", iz=" << src_loc[2] << ")" << std::endl;
}