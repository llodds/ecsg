#include "field.hpp"
#include <iostream>

void compute_res_sxz_vx(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec)
{
    for(int i = 0; i < num_grids-1; i++)
    {
        real res = 0.0;
        for(int i1 = 1; i1 < grids[i].size[1]; i1++)
            for(int i0 = 1; i0 < grids[i].size[0]; i0++)
                res += sim[i][0].GetVal(i0,i1,grids[i].size[2])*
                (sim[i][7].GetVal(i0,i1,grids[i].size[2]-1) + sim[i][7].GetVal(i0,i1,grids[i].size[2]));
        for(int i1 = 1; i1 < grids[i+1].size[1]; i1++)
            for(int i0 = 1; i0 < grids[i+1].size[0]; i0++)
                res -= 4*sim[i+1][0].GetVal(i0,i1,0)*
                (sim[i+1][7].GetVal(i0,i1,-1) + sim[i+1][7].GetVal(i0,i1,0));
        res_vec[i][0].push_back(res);
        res_vec[i][1].push_back(sim[i+1][0].AvgAbs());
        res_vec[i][2].push_back(sim[i+1][7].AvgAbs());
        
        if(fabs(res)>1e-20)
        {
            std::cerr << "==== grid interface " << i << "-" << i+1 << std::endl;
            std::cerr << "\x1B[0m\x1B[31m res of [vx-sxz] = " << res << std::endl;
            std::cerr << " vx2.AvgAbs = " << sim[i+1][0].AvgAbs() << std::endl;
            std::cerr << " sxz2.AvgAbs = " << sim[i+1][7].AvgAbs() << "\x1B[0m" << std::endl;
        }
    }
}


void compute_res_syz_vy(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec)
{
    for(int i = 0; i < num_grids-1; i++)
    {
        real res = 0.0;
        for(int i1 = 0; i1 < grids[i].size[1]; i1++)
            for(int i0 = 0; i0 < grids[i].size[0]; i0++)
                res += sim[i][1].GetVal(i0,i1,grids[i].size[2])*
                (sim[i][6].GetVal(i0,i1,grids[i].size[2]-1) + sim[i][6].GetVal(i0,i1,grids[i].size[2]));
        for(int i1 = 0; i1 < grids[i+1].size[1]; i1++)
            for(int i0 = 0; i0 < grids[i+1].size[0]; i0++)
                res -= 4*sim[i+1][1].GetVal(i0,i1,0)*
                (sim[i+1][6].GetVal(i0,i1,-1) + sim[i+1][6].GetVal(i0,i1,0));
        res_vec[i][3].push_back(res);
        res_vec[i][4].push_back(sim[i+1][1].AvgAbs());
        res_vec[i][5].push_back(sim[i+1][6].AvgAbs());
        
        if(fabs(res)>1e-20)
        {
            std::cerr << "==== grid interface " << i << "-" << i+1 << std::endl;
            std::cerr << "\x1B[0m\x1B[32m res of [vy-syz] = " << res << std::endl;
            std::cerr << " vy2.AvgAbs = " << sim[i+1][1].AvgAbs() << std::endl;
            std::cerr << " syz2.AvgAbs = " << sim[i+1][6].AvgAbs() << "\x1B[0m" << std::endl;
        }
    }
}



void compute_res_szz_vz(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec)
{
    for(int i = 0; i < num_grids-1; i++)
    {
        real res = 0.0;
        for(int i1 = 1; i1 < grids[i].size[1]; i1++)
            for(int i0 = 0; i0 < grids[i].size[0]; i0++)
                res += sim[i][5].GetVal(i0,i1,grids[i].size[2])*
                (sim[i][2].GetVal(i0,i1,grids[i].size[2]-1) + sim[i][2].GetVal(i0,i1,grids[i].size[2]));
        for(int i1 = 1; i1 < grids[i+1].size[1]; i1++)
            for(int i0 = 0; i0 < grids[i+1].size[0]; i0++)
                res -= 4*sim[i+1][5].GetVal(i0,i1,0)*
                (sim[i+1][2].GetVal(i0,i1,-1) + sim[i+1][2].GetVal(i0,i1,0));
        
        res_vec[i][6].push_back(res);
        res_vec[i][7].push_back(sim[i+1][2].AvgAbs());
        res_vec[i][8].push_back(sim[i+1][5].AvgAbs());
        
        if(fabs(res)>1e-20)
        {
            std::cerr << "==== grid interface " << i << "-" << i+1 << std::endl;
            std::cerr << "\x1B[0m\x1B[34m res of [vz-szz] = " << res << std::endl;
            std::cerr << " vz2.AvgAbs = " << sim[i+1][2].AvgAbs() << std::endl;
            std::cerr << " szz2.AvgAbs = " << sim[i+1][5].AvgAbs() << "\x1B[0m" << std::endl;
        }
    }
}