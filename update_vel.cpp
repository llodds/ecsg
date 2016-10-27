#include "field.hpp"
#include <omp.h>

void UpdateVel(Field *sim, const Field &buoy, real **coeff, int radius,
               int top ,int bottom, int top_bc)
{
    Field *vx = sim;
    Field *vy = sim+1;
    Field *vz = sim+2;
    Field *sxx = sim+3;
    Field *syy = sim+4;
    Field *szz = sim+5;
    Field *syz = sim+6;
    Field *sxz = sim+7;
    Field *sxy = sim+8;
    
    if(top && (top_bc == 1)) //free surface boundary condition on the top
        top = 0;
    
    //  Update vx
    //  (Note: shouldn't update vx from vx->GetStart(0) to vx->GetEnd(0)
    //  because txx may not have symmetrix values about vx)
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vx->GetStart(2) + top; iz <= vx->GetEnd(2) - bottom; iz++)
            for(int iy = vx->GetStart(1)+1; iy < vx->GetEnd(1); iy++)
                for(int ir = 0; ir < radius; ir++)
                    for(int ix = vx->GetStart(0)+1; ix < vx->GetEnd(0); ix++)
                        vx->SetVal(ix,iy,iz) += buoy.GetVal(ix,iy,iz)
                        *(coeff[0][ir]*(sxx->GetVal(ix+ir,iy,iz)-sxx->GetVal(ix-ir-1,iy,iz))
                          +coeff[1][ir]*(sxy->GetVal(ix,iy+ir,iz)-sxy->GetVal(ix,iy-ir-1,iz))
                          +coeff[2][ir]*(sxz->GetVal(ix,iy,iz+ir)-sxz->GetVal(ix,iy,iz-ir-1)));
        
        //  Update vy
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vy->GetStart(2) + top; iz <= vy->GetEnd(2) - bottom; iz++)
            for(int iy = vy->GetStart(1); iy <= vy->GetEnd(1); iy++)
                for(int ir = 0; ir < radius; ir++)
                    for(int ix = vy->GetStart(0); ix <= vy->GetEnd(0); ix++)
                        vy->SetVal(ix,iy,iz) += 0.25*(buoy.GetVal(ix,iy,iz)
                                                      +buoy.GetVal(ix+1,iy,iz)
                                                      +buoy.GetVal(ix,iy+1,iz)
                                                      +buoy.GetVal(ix+1,iy+1,iz))
                        *(coeff[0][ir]*(sxy->GetVal(ix+ir+1,iy,iz)-sxy->GetVal(ix-ir,iy,iz))
                          +coeff[2][ir]*(syz->GetVal(ix,iy,iz+ir)-syz->GetVal(ix,iy,iz-ir-1))
                          +coeff[1][ir]*(syy->GetVal(ix,iy+ir+1,iz)-syy->GetVal(ix,iy-ir,iz)));
        
        //  Update vz
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int iy = vz->GetStart(1)+1; iy < vz->GetEnd(1); iy++)
                for(int ir = 0; ir < radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vz->SetVal(ix,iy,iz) += 0.25*(buoy.GetVal(ix,iy,iz)
                                                      +buoy.GetVal(ix,iy,iz+1)
                                                      +buoy.GetVal(ix+1,iy,iz)
                                                      +buoy.GetVal(ix+1,iy,iz+1))
                        *(coeff[0][ir]*(sxz->GetVal(ix+ir+1,iy,iz)-sxz->GetVal(ix-ir,iy,iz))
                          +coeff[2][ir]*(szz->GetVal(ix,iy,iz+ir+1)-szz->GetVal(ix,iy,iz-ir))
                          +coeff[1][ir]*(syz->GetVal(ix,iy+ir,iz)-syz->GetVal(ix,iy-ir-1,iz)));
        
        
        
        //  Update ghost velocities
        //  Case 1: If domain contains z[start]
        if(top)
        {
            int iz = vz->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                for(int ir = 1; ir <= radius - top_bc; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vz->SetVal(ix,iy,iz-ir) = -vz->GetVal(ix,iy,iz+ir-1);
            
            iz = vx->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vx->GetStart(1); iy <= vx->GetEnd(1); iy++)
                for(int ir = 1; ir < radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vx->SetVal(ix,iy,iz-ir) = -vx->GetVal(ix,iy,iz+ir);
            
            iz = vy->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vy->GetStart(1); iy <= vy->GetEnd(1); iy++)
                for(int ir = 1; ir < radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vy->SetVal(ix,iy,iz-ir) = -vy->GetVal(ix,iy,iz+ir);
        }
        
        //  Case 2: If domain contains z[end]
        if(bottom)
        {
            int iz = vz->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                for(int ir = 1; ir <= radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vz->SetVal(ix,iy,iz+ir) = -vz->GetVal(ix,iy,iz-ir+1);
            
            iz = vx->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vx->GetStart(1); iy <= vx->GetEnd(1); iy++)
                for(int ir = 1; ir < radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vx->SetVal(ix,iy,iz+ir) = -vx->GetVal(ix,iy,iz-ir);
            
            iz = vy->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int iy = vy->GetStart(1); iy <= vy->GetEnd(1); iy++)
                for(int ir = 1; ir < radius; ir++)
                    for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                        vy->SetVal(ix,iy,iz+ir) = -vy->GetVal(ix,iy,iz-ir);
        }
        
        //   Case 3: If domain contains x[start]
        int ix = vz->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vz->SetVal(ix-ir,iy,iz) = -vz->GetVal(ix+ir-1,iy,iz);
        
        ix = vx->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vx->GetStart(2); iz <= vx->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vx->SetVal(ix-ir,iy,iz) = -vx->GetVal(ix+ir,iy,iz);
        
        ix = vy->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vy->GetStart(2); iz <= vy->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vy->SetVal(ix-ir,iy,iz) = -vy->GetVal(ix+ir-1,iy,iz);
        
        //   Case 4: If domain contains x[end]
        ix = vz->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vz->SetVal(ix+ir,iy,iz) = -vz->GetVal(ix-ir+1,iy,iz);
        
        ix = vx->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vx->GetStart(2); iz <= vx->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vx->SetVal(ix+ir,iy,iz) = -vx->GetVal(ix-ir,iy,iz);
        
        ix = vy->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vy->GetStart(2); iz <= vy->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
                    vy->SetVal(ix+ir,iy,iz) = -vy->GetVal(ix-ir+1,iy,iz);
        
        //   Case 5: If domain contains y[start]
        int iy = vz->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                    vz->SetVal(ix,iy-ir,iz) = -vz->GetVal(ix,iy+ir,iz);
        
        iy = vx->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int ix = vx->GetStart(0); ix <= vx->GetEnd(0); ix++)
                    vx->SetVal(ix,iy-ir,iz) = -vx->GetVal(ix,iy+ir,iz);
        
        iy = vy->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int ix = vy->GetStart(0); ix <= vy->GetEnd(0); ix++)
                    vy->SetVal(ix,iy-ir,iz) = -vy->GetVal(ix,iy+ir-1,iz);
        
        //   Case 6: If domain contains y[end]
        iy = vz->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
                    vz->SetVal(ix,iy+ir,iz) = -vz->GetVal(ix,iy-ir,iz);
        
        iy = vx->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir < radius; ir++)
                for(int ix = vx->GetStart(0); ix <= vx->GetEnd(0); ix++)
                    vx->SetVal(ix,iy+ir,iz) = -vx->GetVal(ix,iy-ir,iz);
        
        iy = vy->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
            for(int ir = 1; ir <= radius; ir++)
                for(int ix = vy->GetStart(0); ix <= vy->GetEnd(0); ix++)
                    vy->SetVal(ix,iy+ir,iz) = -vy->GetVal(ix,iy-ir+1,iz);
    }
    
}
