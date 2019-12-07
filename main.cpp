/* Based on uniform/composite staggered grid FDTD scheme, this C++ package
 * aims to simulate elastic wave propagations in general anisotropic material.
 *
 * Copyleft 2016 Muhong Zhou.
 *
 * Author: Muhong Zhou < mz10@rice.edu >
 *
 * Naming Convention:
 *   macros are named as: HAPPY, HAPPYENDING
 *   variables are named as: happy, happy_ending
 *   functions are named as: Happy(), HappyEnding()
 *   header files are named as: _HAPPY, _HAPPY_ENDING
 *
 * == Note: the computational domain includes the boundary.  */

#include <omp.h>
#include <iostream>
#include <rsf.hh>
#include <cassert>
#include <vector>
#include <cmath> //fabs()
#include <ctime>
#include "field.hpp"
#include "array.hpp"
#include "sim_utils.hpp"
#include "test_utils.hpp"


int main(int argc, char* argv[])
{
    sf_init(argc, argv);
    oRSF out; //record any output message
    iRSF par(0); //receive any command line input
    
    //check if OpenMP is enabled
    int num_threads;
    par.get("num_threads",num_threads,1);
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#else
    num_threads=1;
#endif
    
    //check if computing the energy
    bool check_energy; par.get("check_energy",check_energy,false);
    
    //read in FD accuracy order (stencil radius)
    int radius; par.get("radius",radius,0); assert(radius>0);
    
    //read in simulation time step
    real dt; par.get("dt",dt,-1.); assert(dt>0);//unit: [ms]
    
    //read in interpolation mode [ECI(1)/II(0)]
    int imode; par.get("imode",imode,1);
    
    //read in verbosity
    bool verb; par.get("verb",verb,false);
    
    //read in boundary condition on the top
    int top_bc; par.get("top_bc", top_bc, 0);
    
    
    
    //=====================================================
    //  INIT GRIDS AND MATERIAL PARAMETER FIELDS
    int num_grids;
    par.get("num_grids",num_grids,1);
    
    Grid* grids; //grids[#grid]
    Field* buoys; //buoys[#grid]
    
    Field **stiffs_obj, **densities_obj, **compliances_obj; //_obj[#grid][#]
    Field ***stiffs_ptr, ***densities_ptr, ***compliances_ptr; //_ptr[#grid][#]
    
    InitParField(num_grids,
                 check_energy,
                 grids, buoys,
                 stiffs_obj, stiffs_ptr,
                 densities_obj, densities_ptr,
                 compliances_obj, compliances_ptr);
    
    
    
    //=====================================================
    //  INIT DYNAMIC(SIMULATION) FIELDS
    //    sim[][0]: vx
    //    sim[][1]: vy
    //    sim[][2]: vz
    //    sim[][3]: sxx
    //    sim[][4]: syy
    //    sim[][5]: szz
    //    sim[][6]: syz
    //    sim[][7]: sxz
    //    sim[][8]: sxy
    Field **sim;
    Field **vel_copy;
    InitSim(num_grids, radius, check_energy, grids,
            sim, vel_copy);
    
    
    
    //=====================================================
    //  INIT TEMPORARY CONTAINER TO COMPUTE THE GHOST DATA
    Array *g_sxz=nullptr, *g_tmp=nullptr, *A_sxz=nullptr, *B_sxz=nullptr;
    Array *A_vz=nullptr, *B_vz=nullptr;
    InitGhostContainer(num_grids, grids, num_threads,
                       g_sxz, g_tmp, A_sxz, B_sxz, A_vz, B_vz);
    
    
    
    //=====================================================
    //  INIT SOURCE
    int src_grid; //source is located on grids[src_grid]
    int src_loc[3]; //source is located on (src_loc[3]) of grids[src_grid]
    std::valarray<float> src_trace; //source data
    InitSrc(src_grid,src_loc,src_trace,
            grids,num_grids, dt);
    int src_len = src_trace.size();
    
    
    
    //=====================================================
    //  COMPUTE FD coefficients
    real*** coeff;
    InitFDcoeff(num_grids, grids, radius, dt, coeff);
    
    
    
    //=====================================================
    //  INIT ENERGY, MOVIES, RECEIVER TRACES, RESIDUE TRACE CONTAINER (verb==true)
    //init energy container
    std::vector<real> energy_vec;
    
    //init movies
    real movie_dt;
    par.get("movie_dt", movie_dt, -1.); assert(movie_dt>0);
    int movie_nt;
    par.get("movie_nt", movie_nt, -1); assert(movie_nt>0);
    bool movie_flags[3];
    par.get("movie_vx_flag", movie_flags[0], false);
    par.get("movie_vy_flag", movie_flags[1], false);
    par.get("movie_vz_flag", movie_flags[2], false);
    
    oRSF*** movie = new oRSF** [3];
    InitMovie(num_grids,movie_flags,sim,grids,movie_nt,movie_dt,movie);
    
    //init receiver traces;
    int n_rec;
    int **rec_loc;
    std::vector<real>* rec_data;
    InitRecTrace(num_grids,grids,n_rec,rec_loc,rec_data);
    
    //init residue traces
    std::vector<float> **res_vec;
    InitRes(num_grids, res_vec);
    
    
    
    
    
    //=====================================================
    clock_t begin=clock();
    clock_t total_time_ghost=0;
    //+-+-+-+-+-+-+-+-+  START SIMULATION
    for(int niter = 0, iframe = 0; niter <= ceil(movie_dt*movie_nt/dt); niter++)
    {
        //monitor how fast an iteration is
        if(verb)
            std::cerr << "niter = " << niter <<  std::endl;
        
        
        //copy vel fields
        if(check_energy)
            for(int i = 0; i < num_grids; i++)
                for(int j = 0; j < 3; j++)
                    vel_copy[i][j] = sim[i][j];
        
        //===update velocities
        for(int i = 0; i < num_grids; i++)
            UpdateVel(sim[i], buoys[i], coeff[i], radius, i==0, i==num_grids-1, top_bc);
        
        //record receiver traces
        for(int i = 0; i < n_rec; i++)
            rec_data[i].push_back(sim[rec_loc[i][0]][0].GetVal(rec_loc[i][1],rec_loc[i][2],rec_loc[i][3]));
        
        //check if (3D.ECI)(a/b) are satisfied --- YES!!!
        if(verb && num_grids>1)
        {
            compute_res_sxz_vx(num_grids, grids, sim, res_vec);
            compute_res_syz_vy(num_grids, grids, sim, res_vec);
        }
        
        //===add source (assume source is only added to vx)
        if(niter < src_len)
            sim[src_grid][0].SetVal(src_loc[0],src_loc[1],src_loc[2]) += 1000*src_trace[niter];
        //===add source (explosive source)
        //if(niter < src_len)
        //        {
        //            sim[src_grid][3].SetVal(src_loc[0],src_loc[1],src_loc[2]) += src_trace[niter];
        //            sim[src_grid][4].SetVal(src_loc[0],src_loc[1],src_loc[2]) += src_trace[niter];
        //            sim[src_grid][5].SetVal(src_loc[0],src_loc[1],src_loc[2]) += src_trace[niter];
        //        }
        
        clock_t begin_ghost=clock();
        //===update ghost vz
        for(int i = 1; i < num_grids; i++)
            UpdateGhostVZ(sim[i-1], sim[i], stiffs_ptr[i-1], stiffs_ptr[i],
                          coeff[i-1], coeff[i],
                          A_vz[i-1], B_vz[i-1], imode);
        total_time_ghost += static_cast<double>(clock()-begin_ghost);
        
        //check if (3D.ECI)(c) is satisfied --- YES!!!
        if(verb && num_grids>1)
            compute_res_szz_vz(num_grids, grids, sim, res_vec);
        
        //compute the energy at this time = niter*dt
        if(check_energy)
        {
            real energy = 0;
            for(int i = 0; i < num_grids; i++)
                energy += ComputeEnergy(sim[i],vel_copy[i],
                                        densities_ptr[i],compliances_ptr[i],
                                        grids[i].step);
            energy_vec.push_back(energy);
        }
        
        //===update stresses
        for(int i = 0; i < num_grids; i++)
            UpdateStress(sim[i], stiffs_ptr[i], coeff[i], radius, i==0, i==num_grids-1, top_bc);
        
        
        
        //===test SBP property
//        if(check_energy)
//            TestSBP(sim);
        
        
        begin_ghost=clock();
        //===update ghost sxz,syz
        for(int i = 1; i < num_grids; i++)
        {
            UpdateGhostSXZ(sim[i-1], sim[i],
                           &buoys[i-1], &buoys[i],
                           coeff[i-1], coeff[i],
                           g_sxz[i-1], A_sxz[i-1], B_sxz[i-1], g_tmp[i-1],
                           imode);
            UpdateGhostSYZ(sim[i-1], sim[i],
                           &buoys[i-1], &buoys[i],
                           coeff[i-1], coeff[i],
                           imode);
        }
        total_time_ghost += static_cast<double>(clock()-begin_ghost);
        
        //output vx to movie
        if(iframe < movie_nt)
            if((niter*dt<=iframe*movie_dt)&&(iframe*movie_dt<niter*dt+dt))
            {
                PrintMovie(num_grids,movie_flags,sim,movie);
                iframe++;
            }
    }
    //+-+-+-+-+-+-+-+-+ END SIMULATION
    double total_time=static_cast<double>(clock()-begin)/CLOCKS_PER_SEC;
    double total_time_g=static_cast<double>(total_time_ghost)/CLOCKS_PER_SEC;
    std::cout << "Total elapsed time " << total_time <<" s." <<std::endl;
    std::cout << "Time spent in computing interface condition: " << total_time_g <<" s." <<std::endl;
    //=====================================================
    
    
    
    
    
    
    //=====================================================
    //  OUTPUT SIMULATION RESULTS
    //output energy
    if(check_energy)
        PrintToRSF(energy_vec, "energy_trace", dt);
    
    //output receiver data trace
    oRSF rec_out("rec_out");
    PrintRecData(n_rec,rec_data,dt,rec_out);
    
    //output 3D.ECI residue
    if(verb && num_grids>1)
    {
        oRSF res_trace("res_trace");
        res_trace.type(SF_FLOAT);
        res_trace.put("n1", static_cast<int>(res_vec[0][0].size()));
        res_trace.put("n2", static_cast<int>(9));
        res_trace.put("n3", num_grids-1);
        
        std::valarray<float> res_trace_vec(res_vec[0][0].size());
        for(int k = 0; k < num_grids-1; k++)
            for(int j = 0; j < 9; j++)
            {
                for(int i = 0; i < res_vec[0][0].size(); i++)
                    res_trace_vec[i] = res_vec[k][j][i];
                res_trace << res_trace_vec;
            }
    }
    
    
    
    //=====================================================
    //  CLEAN UP
    if(verb && num_grids>1)
        DelRes(num_grids, res_vec);
    DelMovie(num_grids,movie_flags,movie);
    DelRecTrace(n_rec,rec_loc,rec_data);
    DelFDCoeff(num_grids, coeff);
    DelGhostContainer(g_sxz,g_tmp,A_sxz,B_sxz,A_vz,B_vz);
    DelSim(num_grids, check_energy,
           sim, vel_copy);
    DelParField(num_grids,
                check_energy,
                grids, buoys,
                stiffs_obj, stiffs_ptr,
                densities_obj, densities_ptr,
                compliances_obj, compliances_ptr);
    
    //=====================================================
    //  THE END
    std::cerr<<"\a\a\a"<<std::endl;
    return 0;
}
