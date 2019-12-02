#include "sim_utils.hpp"
#include <cmath>
#include <cassert>
#include <omp.h>
#include <math.h>

//TODO vectorization
void UpdateGhostVZ(Field *sim1, Field *sim2,
                   Field *stiff1[9], Field *stiff2[9],
                   real **coeff1, real **coeff2,
                   Array &A, Array &B,
                   int mode)
{
    Field& vx1 = sim1[0];
    Field& vy1 = sim1[1];
    Field& vz1 = sim1[2];
    Field& szz1 = sim1[5];
    
    Field& vx2 = sim2[0];
    Field& vy2 = sim2[1];
    Field& vz2 = sim2[2];
    Field& szz2 = sim2[5];
    
    Field *c33_1 = stiff1[2], *c33_2 = stiff2[2];
    Field *c23_1 = stiff1[3], *c23_2 = stiff2[3];
    Field *c13_1 = stiff1[4], *c13_2 = stiff2[4];
    
    int nx2 = vz2.GetEnd(0)+1;
    int ny2 = vz2.GetEnd(1);
    int nz1 = vz1.GetEnd(2)+1;
    
    assert(mode==0 || mode==1);
    //mode =
    //  0: 2nd order accurate intuitive (energy non-conserving) interpolation
    //  1: energy conserving interpolation
    if(mode==0)//intuitive interpolation (2nd order accurate)
    {
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //second order accurate approximation to ghost vz2
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                    vz2.SetVal(i0,i1,-1) = 0.25*(vz1.GetVal(2*i0,2*i1,nz1-2)
                                                 + vz1.GetVal(2*i0+1,2*i1,nz1-2)
                                                 + vz1.GetVal(2*i0,2*i1,nz1-1)
                                                 + vz1.GetVal(2*i0+1,2*i1,nz1-1));
            
            //second order accurate approximation to ghost vz1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2-1; i0++)
                {
                    vz1.SetVal(2*i0+1,2*i1,nz1) = vz1.GetVal(2*i0+1,2*i1,nz1-1)/3
                    + vz2.GetVal(i0,i1,0)/2
                    + vz2.GetVal(i0+1,i1,0)/6;
                    
                    vz1.SetVal(2*i0+2,2*i1,nz1) = vz1.GetVal(2*i0+2,2*i1,nz1-1)/3
                    + vz2.GetVal(i0,i1,0)/6
                    + vz2.GetVal(i0+1,i1,0)/2;
                    
                    vz1.SetVal(2*i0+1,2*i1+1,nz1) = vz1.GetVal(2*i0+1,2*i1+1,nz1-1)/3
                    + vz2.GetVal(i0,i1,0)/4
                    + vz2.GetVal(i0+1,i1,0)/12
                    + vz2.GetVal(i0,i1+1,0)/4
                    + vz2.GetVal(i0+1,i1+1,0)/12;
                    
                    vz1.SetVal(2*i0+2,2*i1+1,nz1) = vz1.GetVal(2*i0+2,2*i1+1,nz1-1)/3
                    + vz2.GetVal(i0,i1,0)/12
                    + vz2.GetVal(i0+1,i1,0)/4
                    + vz2.GetVal(i0,i1+1,0)/12
                    + vz2.GetVal(i0+1,i1+1,0)/4;
                }
            
            //option 1: use extrapolation to update boundary to ensure second order accuracy
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
            {
                int i0 = vz1.GetStart(0);
                vz1.SetVal(i0,2*i1,nz1) = vz1.GetVal(i0,2*i1,nz1-1)*2/3
                - vz1.GetVal(i0+1,2*i1,nz1-1)/3
                + vz2.GetVal(i0,i1,0)*2/3;
                
                vz1.SetVal(i0,2*i1+1,nz1) = vz1.GetVal(i0,2*i1+1,nz1-1)*2/3
                - vz1.GetVal(i0+1,2*i1,nz1-1)/3
                + vz2.GetVal(i0,i1,0)/3
                + vz2.GetVal(i0,i1+1,0)/3;
                
                i0 = vz1.GetEnd(0);
                vz1.SetVal(i0,2*i1,nz1) = vz1.GetVal(i0,2*i1,nz1-1)*2/3
                - vz1.GetVal(i0-1,2*i1,nz1-1)/3
                + vz2.GetVal(i0,i1,0)*2/3;
                
                vz1.SetVal(i0,2*i1+1,nz1) = vz1.GetVal(i0,2*i1+1,nz1-1)*2/3
                - vz1.GetVal(i0-1,2*i1,nz1-1)/3
                + vz2.GetVal(i0,i1,0)/3
                + vz2.GetVal(i0,i1+1,0)/3;
            }
            
            //option 2: do not use extrapolation to update boundary
            //#ifdef _OPENMP
            //#pragma omp for
            //#endif
            //            for(int i1 = 0; i1 < ny2; i1++)
            //                for(int i0 = vz1.GetStart(0); i0 <= vz1.GetEnd(0); i0 += vz1.GetEnd(0)-vz1.GetStart(0))
            //                {
            //                    vz1.SetVal(i0,2*i1,nz1) = vz1.GetVal(i0,2*i1,nz1-1)/3
            //                    + vz2.GetVal(i0,i1,0)*2/3;
            //
            //                    vz1.SetVal(i0,2*i1+1,nz1) = vz1.GetVal(i0,2*i1+1,nz1-1)/3
            //                    + vz2.GetVal(i0,i1,0)/3
            //                    + vz2.GetVal(i0,i1+1,0)/3;
            //                }
        }
    }else if(mode==1)
    {
        //use Tridiagonal Matrix Algorithm (TDMA) to solve
        //each szz2(i+0.5,j,-1)
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int i0 = 0; i0 < nx2; i0++)
        {
#ifdef _OPENMP
            int i_omp = omp_get_thread_num();
#else
            int i_omp = 0;
#endif
            //compute A and b, A*g_tzz2=b
            for(int i1 = 0; i1 < ny2-1; i1++)
            {
                //coeff for szz2(i+0.5,j-1,-1)
                A.SetVal(i1,0,i_omp) = 0.5*coeff1[2][0]*(c33_1->GetVal(2*i0,2*i1+1,nz1)
                                                         +c33_1->GetVal(2*i0+1,2*i1+1,nz1));
                //coeff for szz2(i+0.5,j+1,-1)
                A.SetVal(i1,2,i_omp) = 0.5*coeff1[2][0]*(c33_1->GetVal(2*i0,2*i1+3,nz1)
                                                         +c33_1->GetVal(2*i0+1,2*i1+3,nz1));
                //coeff for szz2(i+0.5,j,-1)
                A.SetVal(i1,1,i_omp) = 8*coeff2[2][0]*c33_2->GetVal(i0,i1+1,0)
                + 2*coeff1[2][0]*(c33_1->GetVal(2*i0,2*i1+2,nz1)+c33_1->GetVal(2*i0+1,2*i1+2,nz1))
                + A.GetVal(i1,0,i_omp) + A.GetVal(i1,2,i_omp);
                
                B.SetVal(i1,i_omp) = 8*(szz2.GetVal(i0,i1+1,0)
                                        + c13_2->GetVal(i0,i1+1,0)*coeff2[0][0]*(vx2.GetVal(i0+1,i1+1,0)
                                                                                 -vx2.GetVal(i0,i1+1,0))
                                        + c23_2->GetVal(i0,i1+1,0)*coeff2[1][0]*(vy2.GetVal(i0,i1+1,0)
                                                                                 -vy2.GetVal(i0,i1,0))
                                        + c33_2->GetVal(i0,i1+1,0)*coeff2[2][0]*vz2.GetVal(i0,i1+1,0))
                - 2*(szz1.GetVal(2*i0,2*i1+2,nz1)
                     + c13_1->GetVal(2*i0,2*i1+2,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+1,2*i1+2,nz1)
                                                                    - vx1.GetVal(2*i0,2*i1+2,nz1))
                     + c23_1->GetVal(2*i0,2*i1+2,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0,2*i1+2,nz1)
                                                                    -vy1.GetVal(2*i0,2*i1+1,nz1))
                     + c33_1->GetVal(2*i0,2*i1+2,nz1)*coeff1[2][0]*(vz2.GetVal(i0,i1+1,0)
                                                                    -2*vz1.GetVal(2*i0,2*i1+2,nz1-1)))
                - 2*(szz1.GetVal(2*i0+1,2*i1+2,nz1)
                     + c13_1->GetVal(2*i0+1,2*i1+2,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+2,2*i1+2,nz1)
                                                                      - vx1.GetVal(2*i0+1,2*i1+2,nz1))
                     + c23_1->GetVal(2*i0+1,2*i1+2,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0+1,2*i1+2,nz1)
                                                                      -vy1.GetVal(2*i0+1,2*i1+1,nz1))
                     + c33_1->GetVal(2*i0+1,2*i1+2,nz1)*coeff1[2][0]*(vz2.GetVal(i0,i1+1,0)
                                                                      -2*vz1.GetVal(2*i0+1,2*i1+2,nz1-1)))
                - (szz1.GetVal(2*i0, 2*i1+1, nz1)
                   + c13_1->GetVal(2*i0,2*i1+1,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+1,2*i1+1,nz1)
                                                                  - vx1.GetVal(2*i0,2*i1+1,nz1))
                   + c23_1->GetVal(2*i0,2*i1+1,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0,2*i1+1,nz1)
                                                                  -vy1.GetVal(2*i0,2*i1,nz1))
                   + c33_1->GetVal(2*i0,2*i1+1,nz1)*coeff1[2][0]*(0.5*(vz2.GetVal(i0,i1+1,0)
                                                                       +vz2.GetVal(i0,i1,0))
                                                                  -2*vz1.GetVal(2*i0,2*i1+1,nz1-1)))
                - (szz1.GetVal(2*i0, 2*i1+3, nz1)
                   + c13_1->GetVal(2*i0,2*i1+3,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+1,2*i1+3,nz1)
                                                                  - vx1.GetVal(2*i0,2*i1+3,nz1))
                   + c23_1->GetVal(2*i0,2*i1+3,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0,2*i1+3,nz1)
                                                                  -vy1.GetVal(2*i0,2*i1+2,nz1))
                   + c33_1->GetVal(2*i0,2*i1+3,nz1)*coeff1[2][0]*(0.5*(vz2.GetVal(i0,i1+2,0)
                                                                       +vz2.GetVal(i0,i1+1,0))
                                                                  -2*vz1.GetVal(2*i0,2*i1+3,nz1-1)))
                - (szz1.GetVal(2*i0+1, 2*i1+1, nz1)
                   + c13_1->GetVal(2*i0+1,2*i1+1,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+2,2*i1+1,nz1)
                                                                    - vx1.GetVal(2*i0+1,2*i1+1,nz1))
                   + c23_1->GetVal(2*i0+1,2*i1+1,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0+1,2*i1+1,nz1)
                                                                    -vy1.GetVal(2*i0+1,2*i1,nz1))
                   + c33_1->GetVal(2*i0+1,2*i1+1,nz1)*coeff1[2][0]*(0.5*(vz2.GetVal(i0,i1+1,0)
                                                                         +vz2.GetVal(i0,i1,0))
                                                                    -2*vz1.GetVal(2*i0+1,2*i1+1,nz1-1)))
                - (szz1.GetVal(2*i0+1,2*i1+3,nz1)
                   + c13_1->GetVal(2*i0+1,2*i1+3,nz1)*coeff1[0][0]*(vx1.GetVal(2*i0+2,2*i1+3,nz1)
                                                                    - vx1.GetVal(2*i0+1,2*i1+3,nz1))
                   + c23_1->GetVal(2*i0+1,2*i1+3,nz1)*coeff1[1][0]*(vy1.GetVal(2*i0+1,2*i1+3,nz1)
                                                                    -vy1.GetVal(2*i0+1,2*i1+2,nz1))
                   + c33_1->GetVal(2*i0+1,2*i1+3,nz1)*coeff1[2][0]*(0.5*(vz2.GetVal(i0,i1+2,0)
                                                                         +vz2.GetVal(i0,i1+1,0))
                                                                    -2*vz1.GetVal(2*i0+1,2*i1+3,nz1-1)));
            }
            
            //=====================================================
            //=====================================================
            //=====================================================
            //compute ghost szz2 using TDMA
            //1. forward elimination
            for(int i1 = 1; i1 < ny2-1; i1++)
            {
                real m = A.GetVal(i1,0,i_omp)/A.GetVal(i1-1,1,i_omp);
                A.SetVal(i1,1,i_omp) -= m*A.GetVal(i1-1,2,i_omp);
                B.SetVal(i1,i_omp) -= m*B.GetVal(i1-1,i_omp);
            }
            //2. backward elimination
            vz2.SetVal(i0,ny2-1,-1) = B.GetVal(ny2-2,i_omp)/A.GetVal(ny2-2,1,i_omp);
            for(int i1 = ny2-2; i1 >= 1; i1--)
                vz2.SetVal(i0,i1,-1) = (B.GetVal(i1-1,i_omp) - A.GetVal(i1-1,2,i_omp)*vz2.GetVal(i0,i1+1,-1))/A.GetVal(i1-1,1,i_omp);
        }
        
        //compute ghost szz1
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                {
                    real tmp = vz2.GetVal(i0,i1,0) + vz2.GetVal(i0,i1,-1);
                    vz1.SetVal(2*i0,2*i1,nz1) = -vz1.GetVal(2*i0,2*i1,nz1-1) + tmp;
                    vz1.SetVal(2*i0+1,2*i1,nz1) = -vz1.GetVal(2*i0+1,2*i1,nz1-1) + tmp;
                }
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                {
                    real tmp = 0.5*(vz2.GetVal(i0,i1,0) + vz2.GetVal(i0,i1,-1)
                                    + vz2.GetVal(i0,i1+1,0) + vz2.GetVal(i0,i1+1,-1));
                    vz1.SetVal(2*i0,2*i1+1,nz1) = -vz1.GetVal(2*i0,2*i1+1,nz1-1) + tmp;
                    vz1.SetVal(2*i0+1,2*i1+1,nz1) = -vz1.GetVal(2*i0+1,2*i1+1,nz1-1) + tmp;
                }
        }
    }
}
