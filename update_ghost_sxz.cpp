#include "sim_utils.hpp"
#include <algorithm>    // std::min
#include <cassert>
#include <omp.h>
#include "array.hpp"

/*
 sim1: ptr to dynamic fields on the fine grid
 buoy1: ptr to buoyancy fields on the fine grid
 sim2: ptr to dynamic fields on the coarse grid
 buoy2: ptr to buoyancy fields on the coarse grid
 
 */

//Assume radius=1, H=2h
void UpdateGhostSXZ(Field *sim1, Field *sim2,
                    Field *buoy1, Field *buoy2,
                    real **coeff1, real **coeff2,
                    Array &g_sxz2, Array &A, Array &B,
                    Array &g_tmp,
                    int mode)
{
    Field *vx1 = sim1,  *vx2 = sim2;
    Field *sxx1 = sim1+3, *sxx2 = sim2+3;
    Field *sxz1 = sim1+7, *sxz2 = sim2+7;
    Field *sxy1 = sim1+8, *sxy2 = sim2+8;
    
    int nx2, ny2; //g_sxz2[0:ny2][0:nx2]
    g_sxz2.size(&nx2, &ny2); nx2--; ny2--;
    int nz1 = buoy1->GetEnd(2);
    
    assert(mode==0 || mode==1);
    //mode =
    //  0: 2nd order accurate intuitive (energy non-conserving) interpolation
    //  1: energy conserving interpolation
    if(mode==0)//intuitive interpolation (2nd order accurate)
    {
        //update ghost sxz2
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //second order accurate approximation to ghost sxz2
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz2->SetVal(i0,i1,-1) = 0.5*(sxz1->GetVal(2*i0,2*i1,nz1-2)
                                                  + sxz1->GetVal(2*i0,2*i1,nz1-1));
            //second order accurate approximation to ghost sxz1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0,2*i1,nz1) = sxz1->GetVal(2*i0,2*i1,nz1-1)/3
                    + sxz2->GetVal(i0,i1,0)*2/3.0;
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0+1,2*i1,nz1) = sxz1->GetVal(2*i0+1,2*i1,nz1-1)/3
                    + sxz2->GetVal(i0,i1,0)/3 + sxz2->GetVal(i0+1,i1,0)/3;
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0,2*i1+1,nz1) = sxz1->GetVal(2*i0,2*i1+1,nz1-1)/3
                    + sxz2->GetVal(i0,i1,0)/3 + sxz2->GetVal(i0,i1+1,0)/3;
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0+1,2*i1+1,nz1) = sxz1->GetVal(2*i0+1,2*i1+1,nz1-1)/3
                    + sxz2->GetVal(i0,i1,0)/6 + sxz2->GetVal(i0,i1+1,0)/6
                    + sxz2->GetVal(i0+1,i1,0)/6 + sxz2->GetVal(i0+1,i1+1,0)/6;
        }
        
    }else if(mode==1)//energy-conserving interpolation
    {
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //compute A
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2-1; i1++)
                for(int i0 = 0; i0 < nx2-1; i0++)
                {
                    A.SetVal(i0,i1,0) = (i1!=0)*(i0!=0)*0.25*buoy1->GetVal(2*i0+1,2*i1+1,nz1)*coeff1[2][0];
                    
                    A.SetVal(i0,i1,1) = (i0!=0)*coeff1[2][0]*(buoy1->GetVal(2*i0+1,2*i1+2,nz1)
                                                              + 0.25*buoy1->GetVal(2*i0+1,2*i1+3,nz1)
                                                              + 0.25*buoy1->GetVal(2*i0+1,2*i1+1,nz1));
                    
                    A.SetVal(i0,i1,2) = (i1!=ny2-2)*(i0!=0)*0.25*buoy1->GetVal(2*i0+1,2*i1+3,nz1)*coeff1[2][0];
                    
                    A.SetVal(i0,i1,3) = (i1!=0)*coeff1[2][0]*(buoy1->GetVal(2*i0+2,2*i1+1,nz1)
                                                              + 0.25*buoy1->GetVal(2*i0+3,2*i1+1,nz1)
                                                              + 0.25*buoy1->GetVal(2*i0+1,2*i1+1,nz1));
                    
                    A.SetVal(i0,i1,4) = 16*buoy2->GetVal(i0+1,i1+1,0)*coeff2[2][0]
                    + coeff1[2][0]*(4*buoy1->GetVal(2*i0+2,2*i1+2,nz1)
                                    + buoy1->GetVal(2*i0+3,2*i1+2,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+2,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1+3,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1+1,nz1)
                                    + 0.25*buoy1->GetVal(2*i0+3,2*i1+3,nz1)
                                    + 0.25*buoy1->GetVal(2*i0+1,2*i1+3,nz1)
                                    + 0.25*buoy1->GetVal(2*i0+3,2*i1+1,nz1)
                                    + 0.25*buoy1->GetVal(2*i0+1,2*i1+1,nz1));
                    
                    A.SetVal(i0,i1,5) = (i1!=ny2-2)*coeff1[2][0]*(buoy1->GetVal(2*i0+2,2*i1+3,nz1)
                                                                  + 0.25*buoy1->GetVal(2*i0+3,2*i1+3,nz1)
                                                                  + 0.25*buoy1->GetVal(2*i0+1,2*i1+3,nz1));
                    
                    A.SetVal(i0,i1,6) = (i1!=0)*(i0!=nx2-2)*0.25*buoy1->GetVal(2*i0+3,2*i1+1,nz1)*coeff1[2][0];
                    
                    A.SetVal(i0,i1,7) = (i0!=nx2-2)*coeff1[2][0]*(buoy1->GetVal(2*i0+3,2*i1+2,nz1)
                                                                  + 0.25*buoy1->GetVal(2*i0+3,2*i1+3,nz1)
                                                                  + 0.25*buoy1->GetVal(2*i0+3,2*i1+1,nz1));
                    
                    A.SetVal(i0,i1,8) = (i1!=ny2-2)*(i0!=nx2-2)*coeff1[2][0]*(0.25*buoy1->GetVal(2*i0+3,2*i1+3,nz1));
                }
            
            //compute B
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2-1; i1++)
                for(int i0 = 0; i0 < nx2-1; i0++)
                    B.SetVal(i0,i1) =
                    16*(vx2->GetVal(i0+1,i1+1,0)
                        + buoy2->GetVal(i0+1,i1+1,0)*(coeff2[0][0]*(sxx2->GetVal(i0+1,i1+1,0)
                                                                    -sxx2->GetVal(i0,i1+1,0))
                                                      +coeff2[1][0]*(sxy2->GetVal(i0+1,i1+1,0)
                                                                     -sxy2->GetVal(i0+1,i1,0))
                                                      +coeff2[2][0]*sxz2->GetVal(i0+1,i1+1,0)))
                    - 4*(vx1->GetVal(2*i0+2,2*i1+2,nz1)
                         + buoy1->GetVal(2*i0+2,2*i1+2,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+2,2*i1+2,nz1)
                                                                           -sxx1->GetVal(2*i0+1,2*i1+2,nz1))
                                                             +coeff1[1][0]*(sxy1->GetVal(2*i0+2,2*i1+2,nz1)
                                                                            -sxy1->GetVal(2*i0+2,2*i1+1,nz1))
                                                             +coeff1[2][0]*(sxz2->GetVal(i0+1,i1+1,0)
                                                                            -2*sxz1->GetVal(2*i0+2,2*i1+2,nz1-1))))
                    - 2*(vx1->GetVal(2*i0+3,2*i1+2,nz1)
                         + buoy1->GetVal(2*i0+3,2*i1+2,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+3,2*i1+2,nz1)
                                                                           -sxx1->GetVal(2*i0+2,2*i1+2,nz1))
                                                             +coeff1[1][0]*(sxy1->GetVal(2*i0+3,2*i1+2,nz1)
                                                                            -sxy1->GetVal(2*i0+3,2*i1+1,nz1))
                                                             +coeff1[2][0]*(0.5*(sxz2->GetVal(i0+1,i1+1,0)
                                                                                 +(i0!=nx2-2)*sxz2->GetVal(i0+2,i1+1,0))
                                                                            -2*sxz1->GetVal(2*i0+3,2*i1+2,nz1-1))))
                    - 2*(vx1->GetVal(2*i0+1,2*i1+2,nz1)
                         + buoy1->GetVal(2*i0+1,2*i1+2,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+1,2*i1+2,nz1)
                                                                           -sxx1->GetVal(2*i0,2*i1+2,nz1))
                                                             +coeff1[1][0]*(sxy1->GetVal(2*i0+1,2*i1+2,nz1)
                                                                            -sxy1->GetVal(2*i0+1,2*i1+1,nz1))
                                                             +coeff1[2][0]*(0.5*((i0!=0)*sxz2->GetVal(i0,i1+1,0)
                                                                                 +sxz2->GetVal(i0+1,i1+1,0))
                                                                            -2*sxz1->GetVal(2*i0+1,2*i1+2,nz1-1))))
                    - 2*(vx1->GetVal(2*i0+2,2*i1+3,nz1)
                         + buoy1->GetVal(2*i0+2,2*i1+3,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+2,2*i1+3,nz1)
                                                                           -sxx1->GetVal(2*i0+1,2*i1+3,nz1))
                                                             +coeff1[1][0]*(sxy1->GetVal(2*i0+2,2*i1+3,nz1)
                                                                            -sxy1->GetVal(2*i0+2,2*i1+2,nz1))
                                                             +coeff1[2][0]*(0.5*(sxz2->GetVal(i0+1,i1+1,0)
                                                                                 +(i1!=ny2-2)*sxz2->GetVal(i0+1,i1+2,0))
                                                                            -2*sxz1->GetVal(2*i0+2,2*i1+3,nz1-1))))
                    - 2*(vx1->GetVal(2*i0+2,2*i1+1,nz1)
                         + buoy1->GetVal(2*i0+2,2*i1+1,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+2,2*i1+1,nz1)
                                                                           -sxx1->GetVal(2*i0+1,2*i1+1,nz1))
                                                             +coeff1[1][0]*(sxy1->GetVal(2*i0+2,2*i1+1,nz1)
                                                                            -sxy1->GetVal(2*i0+2,2*i1,nz1))
                                                             +coeff1[2][0]*(0.5*(sxz2->GetVal(i0+1,i1+1,0)
                                                                                 +(i1!=0)*sxz2->GetVal(i0+1,i1,0))
                                                                            -2*sxz1->GetVal(2*i0+2,2*i1+1,nz1-1))))
                    - (vx1->GetVal(2*i0+3,2*i1+3,nz1)
                       + buoy1->GetVal(2*i0+3,2*i1+3,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+3,2*i1+3,nz1)
                                                                         -sxx1->GetVal(2*i0+2,2*i1+3,nz1))
                                                           +coeff1[1][0]*(sxy1->GetVal(2*i0+3,2*i1+3,nz1)
                                                                          -sxy1->GetVal(2*i0+3,2*i1+2,nz1))
                                                           +coeff1[2][0]*(0.25*(sxz2->GetVal(i0+1,i1+1,0)
                                                                                +(i0!=nx2-2)*sxz2->GetVal(i0+2,i1+1,0)
                                                                                +(i1!=ny2-2)*sxz2->GetVal(i0+1,i1+2,0)
                                                                                +(i1!=ny2-2)*(i0!=nx2-2)*sxz2->GetVal(i0+2,i1+2,0))
                                                                          -2*sxz1->GetVal(2*i0+3,2*i1+3,nz1-1))))
                    - (vx1->GetVal(2*i0+1,2*i1+3,nz1)
                       + buoy1->GetVal(2*i0+1,2*i1+3,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+1,2*i1+3,nz1)
                                                                         -sxx1->GetVal(2*i0,2*i1+3,nz1))
                                                           +coeff1[1][0]*(sxy1->GetVal(2*i0+1,2*i1+3,nz1)
                                                                          -sxy1->GetVal(2*i0+1,2*i1+2,nz1))
                                                           +coeff1[2][0]*(0.25*((i0!=0)*sxz2->GetVal(i0,i1+1,0)
                                                                                +sxz2->GetVal(i0+1,i1+1,0)
                                                                                +(i0!=0)*(i1!=ny2-2)*sxz2->GetVal(i0,i1+2,0)
                                                                                +(i1!=ny2-2)*sxz2->GetVal(i0+1,i1+2,0))
                                                                          -2*sxz1->GetVal(2*i0+1,2*i1+3,nz1-1))))
                    - (vx1->GetVal(2*i0+3,2*i1+1,nz1)
                       + buoy1->GetVal(2*i0+3,2*i1+1,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+3,2*i1+1,nz1)
                                                                         -sxx1->GetVal(2*i0+2,2*i1+1,nz1))
                                                           +coeff1[1][0]*(sxy1->GetVal(2*i0+3,2*i1+1,nz1)
                                                                          -sxy1->GetVal(2*i0+3,2*i1,nz1))
                                                           +coeff1[2][0]*(0.25*((i1!=0)*sxz2->GetVal(i0+1,i1,0)
                                                                                +(i1!=0)*(i0!=nx2-2)*sxz2->GetVal(i0+2,i1,0)
                                                                                +sxz2->GetVal(i0+1,i1+1,0)
                                                                                +(i0!=nx2-2)*sxz2->GetVal(i0+2,i1+1,0))
                                                                          -2*sxz1->GetVal(2*i0+3,2*i1+1,nz1-1))))
                    - (vx1->GetVal(2*i0+1,2*i1+1,nz1)
                       + buoy1->GetVal(2*i0+1,2*i1+1,nz1)*(coeff1[0][0]*(sxx1->GetVal(2*i0+1,2*i1+1,nz1)
                                                                         -sxx1->GetVal(2*i0,2*i1+1,nz1))
                                                           +coeff1[1][0]*(sxy1->GetVal(2*i0+1,2*i1+1,nz1)
                                                                          -sxy1->GetVal(2*i0+1,2*i1,nz1))
                                                           +coeff1[2][0]*(0.25*((i0!=0)*(i1!=0)*sxz2->GetVal(i0,i1,0)
                                                                                +(i1!=0)*sxz2->GetVal(i0+1,i1,0)
                                                                                +(i0!=0)*sxz2->GetVal(i0,i1+1,0)
                                                                                +sxz2->GetVal(i0+1,i1+1,0))
                                                                          -2*sxz1->GetVal(2*i0+1,2*i1+1,nz1-1))));
        }
        
        
        //====== use Jacobi iteration method to solve for ghost sxz2 ======
#ifdef _SPFP
        Jacobi2D(g_sxz2, A, B, nx2-1, ny2-1, std::min(1e-3,B.AvgAbs()*1e-4), g_tmp);
#endif
#ifdef _DPFP
        Jacobi2D(g_sxz2, A, B, nx2-1, ny2-1, std::min(1e-12,B.AvgAbs()*1e-5), g_tmp);
#endif
        //=================================================================
        
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //return ghost sxz2
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz2->SetVal(i0,i1,-1) = g_sxz2.GetVal(i0,i1);
            
            //return ghost sxz1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0,2*i1,nz1) = - sxz1->GetVal(2*i0,2*i1,nz1-1)
                    + sxz2->GetVal(i0,i1,-1) + sxz2->GetVal(i0,i1,0);
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 1; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0+1,2*i1,nz1) = - sxz1->GetVal(2*i0+1,2*i1,nz1-1)
                    + 0.5*((i0!=0)*(sxz2->GetVal(i0,i1,-1)
                                    + sxz2->GetVal(i0,i1,0))
                           + (i0!=nx2-1)*(sxz2->GetVal(i0+1,i1,-1)
                                          + sxz2->GetVal(i0+1,i1,0)));
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 1; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0,2*i1+1,nz1) = -sxz1->GetVal(2*i0,2*i1+1,nz1-1)
                    + 0.5*((i1!=0)*(sxz2->GetVal(i0,i1,-1)
                                    + sxz2->GetVal(i0,i1,0))
                           + (i1!=ny2-1)*(sxz2->GetVal(i0,i1+1,-1)
                                          + sxz2->GetVal(i0,i1+1,0)));
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = 0; i1 < ny2; i1++)
                for(int i0 = 0; i0 < nx2; i0++)
                    sxz1->SetVal(2*i0+1,2*i1+1,nz1) = -sxz1->GetVal(2*i0+1,2*i1+1,nz1-1)
                    + 0.25*((i1!=0)*((i0!=0)*(sxz2->GetVal(i0,i1,-1)
                                              + sxz2->GetVal(i0,i1,0))
                                     + (i0!=nx2-1)*(sxz2->GetVal(i0+1,i1,-1)
                                                    + sxz2->GetVal(i0+1,i1,0)))
                            + (i1!=ny2-1)*((i0!=0)*(sxz2->GetVal(i0,i1+1,-1)
                                                    + sxz2->GetVal(i0,i1+1,0))
                                           + (i0!=nx2-1)*(sxz2->GetVal(i0+1,i1+1,-1)
                                                          + sxz2->GetVal(i0+1,i1+1,0))));
        }
    }
}







/*  Perform jacobi iteration method on A(*)X = B and solve for X (EQ(58) in ESG paper)
 *  Input:
 *      X,X_tmp[0:n1+1][0:n0+1]
 *      A[9][0:n1-1][0:n0-1]
 *      B[0:n1-1][0:n0-1]
 *  Output:
 *      X[0:n1+1][0:n0+1]
 */
void Jacobi2D(Array &X, Array &A, Array &B, int n0, int n1, real tol, Array &X_tmp)
{
    //TODO delete
    //real res=0;
    do{
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1=1; i1<=n1; i1++)
                for(int i0=1; i0<=n0; i0++)
                    X_tmp.SetVal(i0,i1) = X.SetVal(i0,i1);
            
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1=1; i1<=n1; i1++)
                for(int i0=1; i0<=n0; i0++)
                    X.SetVal(i0,i1) = (B.GetVal(i0-1,i1-1)
                                       - A.GetVal(i0-1,i1-1,0)*X_tmp.GetVal(i0-1,i1-1)
                                       - A.GetVal(i0-1,i1-1,1)*X_tmp.GetVal(i0-1,i1)
                                       - A.GetVal(i0-1,i1-1,2)*X_tmp.GetVal(i0-1,i1+1)
                                       - A.GetVal(i0-1,i1-1,3)*X_tmp.GetVal(i0,i1-1)
                                       - A.GetVal(i0-1,i1-1,5)*X_tmp.GetVal(i0,i1+1)
                                       - A.GetVal(i0-1,i1-1,6)*X_tmp.GetVal(i0+1,i1-1)
                                       - A.GetVal(i0-1,i1-1,7)*X_tmp.GetVal(i0+1,i1)
                                       - A.GetVal(i0-1,i1-1,8)*X_tmp.GetVal(i0+1,i1+1))/A.GetVal(i0-1,i1-1,4);
        }
        
//        res = Norm2D(X,X_tmp,n0,n1,2);
//        std::cerr<< "---->res=" << res << std::endl;
//        std::cerr<< "---->tol=" << tol << std::endl;
    }while(Norm2D(X,X_tmp,n0,n1,2)>tol);
}

/* compute relative difference between two matrices
 * X1,X2[0:n1-1][0:n0-1]
 * return: ||x1 - x2||
 * mode: [0]: infinity norm, [2]: two norm
 */
real Norm2D(Array &X1, Array &X2, int n0, int n1, int mode)
{
    assert(mode==0 || mode==2);
    real res = 0;
    
    if(mode==2)
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:res)
#endif
        for(int i1 = 1; i1 <= n1; i1++)
            for(int i0 = 1; i0 <= n0; i0++)
            {
                res += (X1.GetVal(i0,i1) - X2.GetVal(i0,i1))*(X1.GetVal(i0,i1) - X2.GetVal(i0,i1));
            }
        
        res = sqrt(res);
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:res)
#endif
        for(int i1 = 1; i1 <= n1; i1++)
            for(int i0 = 1; i0 <= n0; i0++)
            {
                res = std::max(res, static_cast<real>(fabs(X1.GetVal(i0,i1) - X2.GetVal(i0,i1))));
            }
    }
    
    return res;
}