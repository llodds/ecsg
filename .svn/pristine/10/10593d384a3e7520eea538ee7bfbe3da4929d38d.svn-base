#include "sim_utils.hpp"
#include <cassert>
#include <omp.h>

//assume radius=1, H=2h
void UpdateGhostSYZ(Field *sim1, Field *sim2,
                    Field *buoy1, Field *buoy2,
                    real **coeff1, real **coeff2,
                    int mode)
{
    Field *syy1 = sim1+4, *syy2 = sim2+4;
    Field *syz1 = sim1+6, *syz2 = sim2+6;
    Field *sxy1 = sim1+8, *sxy2 = sim2+8;
    
    Field *vy1 = sim1+1, *vy2 = sim2+1;
    
    int nz1 = buoy1->GetEnd(2);
    
    assert(mode==0 || mode==1);
    //mode =
    //  0: 2nd order accurate intuitive (energy non-conserving) interpolation
    //  1: energy conserving interpolation
    if(mode==0)
    {
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //second order accurate approximation to ghost syz2
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz2->GetStart(1); i1 <= syz2->GetEnd(1); i1++)
                for(int i0 = syz2->GetStart(0); i0 <= syz2->GetEnd(0); i0++)
                    syz2->SetVal(i0,i1,-1) = 0.125*(syz1->GetVal(2*i0,2*i1,nz1-2)
                                                    + syz1->GetVal(2*i0+1,2*i1,nz1-2)
                                                    + syz1->GetVal(2*i0,2*i1+1,nz1-2)
                                                    + syz1->GetVal(2*i0+1,2*i1+1,nz1-2)
                                                    + syz1->GetVal(2*i0,2*i1,nz1-1)
                                                    + syz1->GetVal(2*i0+1,2*i1,nz1-1)
                                                    + syz1->GetVal(2*i0,2*i1+1,nz1-1)
                                                    + syz1->GetVal(2*i0+1,2*i1+1,nz1-1));
            
            
            //second order accurate approximation to ghost syz1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz2->GetStart(1); i1 < syz2->GetEnd(1); i1++)
                for(int i0 = syz2->GetStart(0); i0 < syz2->GetEnd(0); i0++)
                {
                    syz1->SetVal(2*i0+1,2*i1+1,nz1) = syz1->GetVal(2*i0+1,2*i1+1,nz1-1)/3
                    + 0.375*syz2->GetVal(i0,i1,0)
                    + 0.125*syz2->GetVal(i0+1,i1,0)
                    + 0.125*syz2->GetVal(i0,i1+1,0)
                    + syz2->GetVal(i0+1,i1+1,0)/24;
                    
                    syz1->SetVal(2*i0+2,2*i1+1,nz1) = syz1->GetVal(2*i0+2,2*i1+1,nz1-1)/3
                    + 0.125*syz2->GetVal(i0,i1,0)
                    + 0.375*syz2->GetVal(i0+1,i1,0)
                    + syz2->GetVal(i0,i1+1,0)/24
                    + 0.125*syz2->GetVal(i0+1,i1+1,0);
                    
                    syz1->SetVal(2*i0+1,2*i1+2,nz1) = syz1->GetVal(2*i0+1,2*i1+2,nz1-1)/3
                    + 0.125*syz2->GetVal(i0,i1,0)
                    + syz2->GetVal(i0+1,i1,0)/24
                    + 0.375*syz2->GetVal(i0,i1+1,0)
                    + 0.125*syz2->GetVal(i0+1,i1+1,0);
                    
                    syz1->SetVal(2*i0+2,2*i1+2,nz1) = syz1->GetVal(2*i0+2,2*i1+2,nz1-1)/3
                    + syz2->GetVal(i0,i1,0)/24
                    + 0.125*syz2->GetVal(i0+1,i1,0)
                    + 0.125*syz2->GetVal(i0,i1+1,0)
                    + 0.375*syz2->GetVal(i0+1,i1+1,0);
                }
            
            //  option1: use extrapolation to update boundary points to ensure second order accuracy
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i0 = syz1->GetStart(0); i0 <= syz1->GetEnd(0); i0++) //update top and bottom edges
            {
                if(i0 == syz1->GetStart(0))
                {
                    syz1->SetVal(i0,syz1->GetStart(1),nz1) = syz1->GetVal(i0,syz1->GetStart(1),nz1-1)
                    - syz1->GetVal(i0+1,syz1->GetStart(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetStart(1)+1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetStart(1)/2,0)*2/3;
                    
                    syz1->SetVal(i0,syz1->GetEnd(1),nz1) = syz1->GetVal(i0,syz1->GetEnd(1),nz1-1)
                    - syz1->GetVal(i0+1,syz1->GetEnd(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetEnd(1)-1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetEnd(1)/2,0)*2/3;
                }else if(i0 == syz1->GetEnd(0))
                {
                    syz1->SetVal(i0,syz1->GetStart(1),nz1) = syz1->GetVal(i0,syz1->GetStart(1),nz1-1)
                    - syz1->GetVal(i0-1,syz1->GetStart(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetStart(1)+1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetStart(1)/2,0)*2/3;
                    
                    syz1->SetVal(i0,syz1->GetEnd(1),nz1) = syz1->GetVal(i0,syz1->GetEnd(1),nz1-1)
                    - syz1->GetVal(i0-1,syz1->GetEnd(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetEnd(1)-1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetEnd(1)/2,0)*2/3;
                }else
                {
                    syz1->SetVal(i0,syz1->GetStart(1),nz1) = syz1->GetVal(i0,syz1->GetStart(1),nz1-1)/3
                    + syz1->GetVal(i0-1+2*(i0%2),syz1->GetStart(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetStart(1)+1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetStart(1)/2,0)*2/3;
                    
                    syz1->SetVal(i0,syz1->GetEnd(1),nz1) = syz1->GetVal(i0,syz1->GetEnd(1),nz1-1)/3
                    + syz1->GetVal(i0-1+2*(i0%2),syz1->GetEnd(1),nz1-1)/3
                    - syz1->GetVal(i0,syz1->GetEnd(1)-1,nz1-1)/3
                    + syz2->GetVal(i0/2,syz1->GetEnd(1)/2,0)*2/3;
                }
            }
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz1->GetStart(1)+1; i1 < syz1->GetEnd(1); i1++) //update left and right edges
            {
                syz1->SetVal(syz1->GetStart(0),i1,nz1) = syz1->GetVal(syz1->GetStart(0),i1,nz1-1)/3
                + syz1->GetVal(syz1->GetStart(0),i1-1+2*(i1%2),nz1-1)/3
                - syz1->GetVal(syz1->GetStart(0)+1,i1,nz1-1)/3
                + syz2->GetVal(syz1->GetStart(0)/2,i1/2,0)*2/3;
                
                syz1->SetVal(syz1->GetEnd(0),i1,nz1) = syz1->GetVal(syz1->GetEnd(0),i1,nz1-1)/3
                + syz1->GetVal(syz1->GetEnd(0),i1-1+2*(i1%2),nz1-1)/3
                - syz1->GetVal(syz1->GetEnd(0)-1,i1,nz1-1)/3
                + syz2->GetVal(syz1->GetEnd(0)/2,i1/2,0)*2/3;
            }
            
            
            //  option 2: do not use extrapolation near the boundary
            //#ifdef _OPENMP
            //#pragma omp for
            //#endif
            //                for(int i0 = syz1->GetStart(0); i0 <= syz1->GetEnd(0); i0++)
            //                {
            //                syz1->SetVal(i0,syz1->GetStart(1),nz1) = syz1->GetVal(i0,syz1->GetStart(1),nz1-1)/3 + 2*syz2->GetVal(i0/2,syz1->GetStart(1)/2,0)/3;
            //                syz1->SetVal(i0,syz1->GetEnd(1),nz1) = syz1->GetVal(i0,syz1->GetEnd(1),nz1-1)/3 + 2*syz2->GetVal(i0/2,syz1->GetEnd(1)/2,0)/3;
            //                 }
            //#ifdef _OPENMP
            //#pragma omp for
            //#endif
            //            for(int i1 = syz1->GetStart(1)+1; i1 < syz1->GetEnd(1); i1++)
            //                for(int i0 = syz1->GetStart(0); i0 <= syz1->GetEnd(0); i0 += syz1->GetEnd(0)-syz1->GetStart(0))
            //                    syz1->SetVal(i0,i1,nz1) = syz1->GetVal(i0,i1,nz1-1)/3 + 2*syz2->GetVal(i0/2,i1/2,0)/3;
        }
    }else if(mode==1)
    {
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
        {
            //TODO vectorize!!!
            //compute ghost syz2
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz2->GetStart(1); i1 <= syz2->GetEnd(1); i1++)
                for(int i0 = syz2->GetStart(0); i0 <= syz2->GetEnd(0); i0++)
                {
                    real c0 = 0.25*(buoy2->GetVal(i0,i1,0)
                                    + buoy2->GetVal(i0+1,i1,0)
                                    + buoy2->GetVal(i0,i1+1,0)
                                    + buoy2->GetVal(i0+1,i1+1,0));
                    real c1 = 0.25*(buoy1->GetVal(2*i0,2*i1,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1,nz1)
                                    + buoy1->GetVal(2*i0,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+1,nz1));
                    real c2 = 0.25*(buoy1->GetVal(2*i0+1,2*i1,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1+1,nz1));
                    real c3 = 0.25*(buoy1->GetVal(2*i0,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0,2*i1+2,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+2,nz1));
                    real c4 = 0.25*(buoy1->GetVal(2*i0+1,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1+1,nz1)
                                    + buoy1->GetVal(2*i0+1,2*i1+2,nz1)
                                    + buoy1->GetVal(2*i0+2,2*i1+2,nz1));
                    
                    syz2->SetVal(i0,i1,-1) = 4*(vy2->GetVal(i0,i1,0)
                                                + c0*(coeff2[0][0]*(sxy2->GetVal(i0+1,i1,0)
                                                                    -sxy2->GetVal(i0,i1,0))
                                                      + coeff2[1][0]*(syy2->GetVal(i0,i1+1,0)
                                                                      -syy2->GetVal(i0,i1,0))
                                                      + coeff2[2][0]*syz2->GetVal(i0,i1,0)))
                    - (vy1->GetVal(2*i0,2*i1,nz1)
                       + c1*(coeff1[0][0]*(sxy1->GetVal(2*i0+1,2*i1,nz1)
                                           -sxy1->GetVal(2*i0,2*i1,nz1))
                             + coeff1[1][0]*(syy1->GetVal(2*i0,2*i1+1,nz1)
                                             -syy1->GetVal(2*i0,2*i1,nz1))
                             + coeff1[2][0]*(syz2->GetVal(i0,i1,0)
                                             - 2*syz1->GetVal(2*i0,2*i1,nz1-1))))
                    - (vy1->GetVal(2*i0+1,2*i1,nz1)
                       + c2*(coeff1[0][0]*(sxy1->GetVal(2*i0+2,2*i1,nz1)
                                           - sxy1->GetVal(2*i0+1,2*i1,nz1))
                             + coeff1[1][0]*(syy1->GetVal(2*i0+1,2*i1+1,nz1)
                                             - syy1->GetVal(2*i0+1,2*i1,nz1))
                             + coeff1[2][0]*(syz2->GetVal(i0,i1,0)
                                             - 2*syz1->GetVal(2*i0+1,2*i1,nz1-1))))
                    - (vy1->GetVal(2*i0,2*i1+1,nz1)
                       + c3*(coeff1[0][0]*(sxy1->GetVal(2*i0+1,2*i1+1,nz1)
                                           - sxy1->GetVal(2*i0,2*i1+1,nz1))
                             + coeff1[1][0]*(syy1->GetVal(2*i0,2*i1+2,nz1)
                                             - syy1->GetVal(2*i0,2*i1+1,nz1))
                             + coeff1[2][0]*(syz2->GetVal(i0,i1,0)
                                             - 2*syz1->GetVal(2*i0,2*i1+1,nz1-1))))
                    - (vy1->GetVal(2*i0+1,2*i1+1,nz1)
                       + c4*(coeff1[0][0]*(sxy1->GetVal(2*i0+2,2*i1+1,nz1)
                                           - sxy1->GetVal(2*i0+1,2*i1+1,nz1))
                             + coeff1[1][0]*(syy1->GetVal(2*i0+1,2*i1+2,nz1)
                                             - syy1->GetVal(2*i0+1,2*i1+1,nz1))
                             + coeff1[2][0]*(syz2->GetVal(i0,i1,0)
                                             - 2*syz1->GetVal(2*i0+1,2*i1+1,nz1-1))));
                    
                    syz2->SetVal(i0,i1,-1) /= (4*c0*coeff2[2][0] + c1*coeff1[2][0]
                                               + c2*coeff1[2][0] + c3*coeff1[2][0]
                                               + c4*coeff1[2][0]);
                }
            
            //compute ghost syz1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz2->GetStart(1); i1 <= syz2->GetEnd(1); i1++)
                for(int i0 = syz2->GetStart(0); i0 <= syz2->GetEnd(0); i0++)
                {
                    real tmp = syz2->GetVal(i0,i1,-1) + syz2->GetVal(i0,i1,0);
                    syz1->SetVal(2*i0,2*i1,nz1) = -syz1->GetVal(2*i0,2*i1,nz1-1) + tmp;
                    syz1->SetVal(2*i0+1,2*i1,nz1) = -syz1->GetVal(2*i0+1,2*i1,nz1-1) + tmp;
                    syz1->SetVal(2*i0,2*i1+1,nz1) = -syz1->GetVal(2*i0,2*i1+1,nz1-1) + tmp;
                    syz1->SetVal(2*i0+1,2*i1+1,nz1) = -syz1->GetVal(2*i0+1,2*i1+1,nz1-1) + tmp;
                }
        }
    }
    
}