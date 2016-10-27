#include <rsf.hh>
#include <vector>
#include "grid.hpp"
#include <cassert>
#include <iostream>


void InitRecTrace(//input
                  int num_grids,
                  Grid grids[],
                  //output
                  int& n_rec,
                  int**& rec_loc,
                  std::vector<real>*& rec_data)
{
    iRSF rec_in("rec_in");
    
    //read in number of receivers
    rec_in.get("n2",n_rec,0);
    
    //read in receiver locations and init receiver grid locations
    rec_loc = new int* [n_rec];
    std::valarray<float> loc_trace(3);
    for(int i = 0; i < n_rec; i++)
    {
        rec_loc[i] = new int [4];
        //rec_loc[0]: #grid
        //rec_loc[1]: ix
        //rec_loc[2]: iy
        //rec_loc[3]: iz
        
        rec_in >> loc_trace;
        real tmp = loc_trace[2];
        rec_loc[i][0] = 0;
        while(tmp > grids[rec_loc[i][0]].step[2]*grids[rec_loc[i][0]].size[2])
        {
            tmp -= grids[rec_loc[i][0]].step[2]*grids[rec_loc[i][0]].size[2];
            rec_loc[i][0]++;
        }
        assert(rec_loc[i][0] < num_grids);
        
        for(int j = 0; j < 2; j++)
            rec_loc[i][j+1] = round((loc_trace[j]-grids[rec_loc[i][0]].origin[j])/grids[rec_loc[i][0]].step[j]);
        rec_loc[i][3] = round(tmp/grids[rec_loc[i][0]].step[2]);
    }
    
    //print out rec_location
    std::cout << "== Record traces at... " << std::endl;
    for(int i = 0; i < n_rec; i++)
        std::cout << "(grid#" << rec_loc[i][0] << ", ix=" << rec_loc[i][1] << ",  iy=" << rec_loc[i][2] << ", iz=" << rec_loc[i][3] << ")" << std::endl;
    
    
    
    //init receiver trace data containers
    rec_data = new std::vector<real> [n_rec];
}


void DelRecTrace(int n_rec,
                 int** rec_loc,
                 std::vector<real>* rec_data)
{
    for(int i = 0; i < n_rec; i++)
        delete [] rec_loc[i];
    delete [] rec_loc;
    
    delete [] rec_data;
}