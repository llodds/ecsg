#include <vector>
#include "field.hpp"
#include "array.hpp"

#ifndef SIM_UTILS
#define SIM_UTILS

//=====================================================
//functions to init and del parameter fields [init_del_par_field.hpp]
void InvertParameters(Field *density_ptr[3], Field density_obj[3],
                      Field *compliance_ptr[9], Field compliance_obj[9],
                      const Field &buoy, Field *stiffness_ptr[9], int aniso_mode);
void InitParField(//input
                  int num_grids,
                  bool check_energy,
                  //output
                  Grid *&grids,
                  Field *&buoys,
                  Field **&stiffs_obj,
                  Field ***&stiffs_ptr,
                  Field **&densities_obj,
                  Field ***&densities_ptr,
                  Field **&compliances_obj,
                  Field ***&compliances_ptr);
void DelParField(//input
                 int num_grids,
                 bool check_energy,
                 //output
                 Grid *&grids,
                 Field *&buoys,
                 Field **&stiffs_obj,
                 Field ***&stiffs_ptr,
                 Field **&densities_obj,
                 Field ***&densities_ptr,
                 Field **&compliances_obj,
                 Field ***&compliances_ptr);


//=====================================================
//functions to init and del simulation fields [init_del_sim.cpp]
void InitSim(int num_grids, int radius, bool check_energy, Grid grids[],
             Field **&sim, Field **&vel_copy);
void DelSim(int num_grids, bool check_energy,
            Field **&sim, Field **&vel_copy);



//=====================================================
//functions to init and del temporary containers used to compute the ghost data [init_del_ghost_container.cpp]
void InitGhostContainer(int num_grids,
                        Grid grids[],
                        int num_threads,
                        Array *& g_sxz,
                        Array *& g_tmp,
                        Array *& A_sxz,
                        Array *& B_sxz,
                        Array *& A_vz,
                        Array *& B_vz);
void DelGhostContainer(Array *& g_sxz,
                       Array *& g_tmp,
                       Array *& A_sxz,
                       Array *& B_sxz,
                       Array *& A_vz,
                       Array *& B_vz);


//=====================================================
//functions to init and del finite difference coefficients fields [init_del_fdcoeff.cpp]
std::vector<real> FDCoeff(int radius);
void InitFDcoeff(int num_grids, Grid grids[], int radius, real dt, real ***& coeff);
void DelFDCoeff(int num_grids, real ***& coeff);


//=====================================================
//function to init source [init_src.cpp]
void InitSrc(//output
             int& src_grid,
             int src_loc[],
             std::valarray<float>& src_trace,
             //input
             Grid grids[],
             int num_grids,
             real dt);


//=====================================================
//function to init and del receiver traces [init_del_rec_trace.cpp]
void InitRecTrace(//input
                  int num_grids,
                  Grid grids[],
                  //output
                  int& n_rec,
                  int**& rec_loc,
                  std::vector<real>*& rec_data);
void DelRecTrace(int n_rec,
                 int** rec_loc,
                 std::vector<real>* rec_data);

//=====================================================
//function to init and del movies [init_del_movie.cpp]
void InitMovie(int num_grids,
               bool movie_flags[3],
               Field **sim,
               Grid grids[],
               int movie_nt,
               real movie_dt,
               oRSF*** &movie);
void InitRSFMovie(const Field &f, int nframe, Grid grid, real time_interval, oRSF *rsf_file);
void DelMovie(int num_grids,
              bool movie_flags[3],
              oRSF*** movie);
void PrintMovie(int num_grids,
                bool movie_flags[3],
                Field **sim,
                oRSF ***movie);

//=====================================================
//function to init and del residue traces [init_del_res_trace.cpp]
void InitRes(int num_grids,
             std::vector<float> **&res_vec);
void DelRes(int num_grids,
            std::vector<float> **&res_vec);



//=====================================================
//functions to update vels/stresses/ghost data [update_*.cpp]
void UpdateVel(Field sim[9], const Field &buoy, real **coeff, int radius, int top ,int bottom, int top_bc);
void UpdateStress(Field sim[9], Field *stiffness_ptr[9], real **coeffs, int radius, int top, int bottom, int top_bc);
void Jacobi2D(Array &X, Array &A, Array &B, int n0, int n1, real tol, Array &X_tmp);
real Norm2D(Array &X1, Array &X2, int n0, int n1, int mode);
void UpdateGhostSXZ(Field sim1[9], Field sim2[9],
                    Field *buoy1, Field *buoy2,
                    real **coeff1, real **coeff2,
                    Array &g_sxz2, Array &A, Array &B,
                    Array &g_tmp, int mode);
void UpdateGhostSYZ(Field sim1[9], Field sim2[9],
                    Field *buoy1, Field *buoy2,
                    real **coeff1, real **coeff2, int mode);
void UpdateGhostVZ(Field sim1[9], Field sim2[9],
                   Field *stiff1[9], Field *stiff2[9],
                   real **coeff1, real **coeff2,
                   Array &A, Array &B, int mode);



//=====================================================
//functions to compute energy, residue traces [compute_*.cpp]
real ComputeEnergy(const Field sim[9], const Field vel_copy[3], Field *density_ptr[3], Field *compliance_ptr[9], const real grid_size[3]);
void compute_res_sxz_vx(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec);
void compute_res_syz_vy(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec);
void compute_res_szz_vz(int num_grids,
                        Grid grids[],
                        Field **sim,
                        std::vector<float> **res_vec);


//=====================================================
//functions to print fields to rsf files. [print.cpp]
void PrintToRSF(const std::vector<real> & val_array, std::string val_name, real dt);
void PrintToRSF(const Field *movie, int nframe,
                Grid grid, real time_interval,
                std::string movie_name);
void PrintFieldToRSF(const Field &f, oRSF *out);
void PrintRecData(int n_rec,
                  std::vector<real>* rec_data,
                  real dt,
                  oRSF &rec_out);

#endif