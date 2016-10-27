#include "array.hpp"
#include "grid.hpp"

void InitGhostContainer(int num_grids,
                        Grid grids[],
                        int num_threads,
                        Array *& g_sxz,
                        Array *& g_tmp,
                        Array *& A_sxz,
                        Array *& B_sxz,
                        Array *& A_vz,
                        Array *& B_vz)
{
    if(num_grids > 1)
    {
        g_sxz = new Array [num_grids-1];
        g_tmp = new Array [num_grids-1];
        A_sxz = new Array [num_grids-1];
        B_sxz = new Array [num_grids-1];
        
        for(int i = 1; i < num_grids; i++)
        {
            g_sxz[i-1] = Array(grids[i].size[0]+1, grids[i].size[1]+1);
            g_tmp[i-1] = Array(grids[i].size[0]+1, grids[i].size[1]+1);
            A_sxz[i-1] = Array(grids[i].size[0]-1, grids[i].size[1]-1,9);
            B_sxz[i-1] = Array(grids[i].size[0]-1, grids[i].size[1]-1);
        }
    }
    
    if(num_grids > 1)
    {
        A_vz = new Array [num_grids-1];
        B_vz = new Array [num_grids-1];
        for(int i = 1; i < num_grids; i++)
        {
            A_vz[i-1] = Array(grids[i].size[1]-1, 3, num_threads);
            B_vz[i-1] = Array(grids[i].size[1]-1, num_threads);
            //note: num_threads is added for OMP parallelization
        }
    }
}


void DelGhostContainer(Array *& g_sxz,
                       Array *& g_tmp,
                       Array *& A_sxz,
                       Array *& B_sxz,
                       Array *& A_vz,
                       Array *& B_vz)
{
    delete [] g_sxz;
    delete [] g_tmp;
    delete [] A_sxz;
    delete [] B_sxz;
    delete [] A_vz;
    delete [] B_vz;
}