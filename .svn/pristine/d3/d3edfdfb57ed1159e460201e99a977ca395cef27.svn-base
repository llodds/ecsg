#include "field.hpp"
#include <omp.h>

void UpdateStress(Field *sim, Field *stiffness_ptr[9], real **coeff, int radius,
                  int top, int bottom, int top_bc)
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
    
    const Field *c11 = stiffness_ptr[0];
    const Field *c22 = stiffness_ptr[1];
    const Field *c33 = stiffness_ptr[2];
    const Field *c23 = stiffness_ptr[3];
    const Field *c13 = stiffness_ptr[4];
    const Field *c12 = stiffness_ptr[5];
    const Field *c44 = stiffness_ptr[6];
    const Field *c55 = stiffness_ptr[7];
    const Field *c66 = stiffness_ptr[8];
    
    
    int flag = 0;
    if(top && top_bc)
        flag = 1;

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
        //  Update szz,sxx,syy
        if(flag){
            int i2 = szz->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = szz->GetStart(1); i1 <= szz->GetEnd(1); i1++)
                for(int ir = 0; ir < radius; ir++)
                    for(int i0 = szz->GetStart(0) ; i0 <= szz->GetEnd(0); i0++)
                    {
                        //  Update sxx
                        sxx->SetVal(i0,i1,i2) +=
                        (c11->GetVal(i0,i1,i2)
                         -c13->GetVal(i0,i1,i2)*c13->GetVal(i0,i1,i2)/c33->GetVal(i0,i1,i2))
                        *coeff[0][ir]*(vx->GetVal(i0+ir+1,i1,i2)-vx->GetVal(i0-ir,i1,i2))
                        +(c12->GetVal(i0,i1,i2)
                          -c13->GetVal(i0,i1,i2)*c23->GetVal(i0,i1,i2)/c33->GetVal(i0,i1,i2))
                        *coeff[1][ir]*(vy->GetVal(i0,i1+ir,i2)-vy->GetVal(i0,i1-ir-1,i2));
                        //  Update syy
                        syy->SetVal(i0,i1,i2) +=
                        (c12->GetVal(i0,i1,i2)
                         -c23->GetVal(i0,i1,i2)*c13->GetVal(i0,i1,i2)/c33->GetVal(i0,i1,i2))
                        *coeff[0][ir]*(vx->GetVal(i0+ir+1,i1,i2)-vx->GetVal(i0-ir,i1,i2))
                        +(c22->GetVal(i0,i1,i2)
                          -c23->GetVal(i0,i1,i2)*c23->GetVal(i0,i1,i2)/c33->GetVal(i0,i1,i2))
                        *coeff[1][ir]*(vy->GetVal(i0,i1+ir,i2)-vy->GetVal(i0,i1-ir-1,i2));
                        //  Update szz (TODO: delete)
                        //szz->SetVal(i0,i1,i2) = 0;
                    }
        }
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = szz->GetStart(2) + flag; i2 <= szz->GetEnd(2); i2++)
            for(int i1 = szz->GetStart(1); i1 <= szz->GetEnd(1); i1++)
                for(int ir = 0; ir < radius; ir++)
                    for(int i0 = szz->GetStart(0); i0 <= szz->GetEnd(0); i0++){
                        //  Update sxx
                        sxx->SetVal(i0,i1,i2) +=
                        c11->GetVal(i0,i1,i2)*coeff[0][ir]
                        *(vx->GetVal(i0+ir+1,i1,i2)-vx->GetVal(i0-ir,i1,i2))
                        +c12->GetVal(i0,i1,i2)*coeff[1][ir]
                        *(vy->GetVal(i0,i1+ir,i2)-vy->GetVal(i0,i1-ir-1,i2))
                        +c13->GetVal(i0,i1,i2)*coeff[2][ir]
                        *(vz->GetVal(i0,i1,i2+ir)-vz->GetVal(i0,i1,i2-ir-1));
                        //  Update syy
                        syy->SetVal(i0,i1,i2) +=
                        c12->GetVal(i0,i1,i2)*coeff[0][ir]
                        *(vx->GetVal(i0+ir+1,i1,i2)-vx->GetVal(i0-ir,i1,i2))
                        +c22->GetVal(i0,i1,i2)*coeff[1][ir]
                        *(vy->GetVal(i0,i1+ir,i2)-vy->GetVal(i0,i1-ir-1,i2))
                        +c23->GetVal(i0,i1,i2)*coeff[2][ir]
                        *(vz->GetVal(i0,i1,i2+ir)-vz->GetVal(i0,i1,i2-ir-1));
                        //  Update szz
                        szz->SetVal(i0,i1,i2) +=
                        c13->GetVal(i0,i1,i2)*coeff[0][ir]
                        *(vx->GetVal(i0+ir+1,i1,i2)-vx->GetVal(i0-ir,i1,i2))
                        +c23->GetVal(i0,i1,i2)*coeff[1][ir]
                        *(vy->GetVal(i0,i1+ir,i2)-vy->GetVal(i0,i1-ir-1,i2))
                        +c33->GetVal(i0,i1,i2)*coeff[2][ir]
                        *(vz->GetVal(i0,i1,i2+ir)-vz->GetVal(i0,i1,i2-ir-1));
                    }
        
        //  Update syz
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = syz->GetStart(2); i2 <= syz->GetEnd(2); i2++)
            for(int i1 = syz->GetStart(1); i1 <= syz->GetEnd(1); i1++)
                for(int ir = 0; ir < radius; ir++)
                    for(int i0 = syz->GetStart(0); i0 <= syz->GetEnd(0); i0++)
                        syz->SetVal(i0,i1,i2) += 0.5*c44->GetVal(i0,i1,i2)
                        *(coeff[2][ir]*(vy->GetVal(i0,i1,i2+ir+1)-vy->GetVal(i0,i1,i2-ir))
                          +coeff[1][ir]*(vz->GetVal(i0,i1+ir+1,i2)-vz->GetVal(i0,i1-ir,i2)));
        
        //  Update sxz
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxz->GetStart(2); i2 <= sxz->GetEnd(2); i2++)
            for(int i1 = sxz->GetStart(1); i1 <= sxz->GetEnd(1); i1++)
                for(int ir = 0; ir < radius; ir++)
                    for(int i0 = sxz->GetStart(0); i0 <= sxz->GetEnd(0); i0++)
                        sxz->SetVal(i0,i1,i2) += 0.5*c55->GetVal(i0,i1,i2)
                        *(coeff[2][ir]*(vx->GetVal(i0,i1,i2+ir+1)-vx->GetVal(i0,i1,i2-ir))
                          +coeff[0][ir]*(vz->GetVal(i0+ir,i1,i2)-vz->GetVal(i0-ir-1,i1,i2)));
        
        //  Update sxy
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxy->GetStart(2); i2 <= sxy->GetEnd(2); i2++)
            for(int i1 = sxy->GetStart(1); i1 <= sxy->GetEnd(1); i1++)
                for(int ir = 0; ir < radius; ir++)
                    for(int i0 = sxy->GetStart(0); i0 <= sxy->GetEnd(0); i0++)
                        sxy->SetVal(i0,i1,i2) += 0.5*c66->GetVal(i0,i1,i2)
                        *(coeff[1][ir]*(vx->GetVal(i0,i1+ir+1,i2)-vx->GetVal(i0,i1-ir,i2))
                          +coeff[0][ir]*(vy->GetVal(i0+ir,i1,i2)-vy->GetVal(i0-ir-1,i1,i2)));
        
        //  Update ghost stresses
        //  If domain contains x[start]
        int i0 = szz->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = szz->GetStart(2); i2 <= szz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = szz->GetStart(1); i1 <= szz->GetEnd(1); i1++)
                    sxx->SetVal(i0-ir,i1,i2) = sxx->GetVal(i0+ir-1,i1,i2);
        
        i0 = sxz->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxz->GetStart(2); i2 <= sxz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = sxz->GetStart(1); i1 <= sxz->GetEnd(1); i1++)
                    sxz->SetVal(i0-ir,i1,i2) = sxz->GetVal(i0+ir,i1,i2);
        
        i0 = sxy->GetStart(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxy->GetStart(2); i2 <= sxy->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = sxy->GetStart(1); i1 <= sxy->GetEnd(1); i1++)
                    sxy->SetVal(i0-ir,i1,i2) = sxy->GetVal(i0+ir,i1,i2);
        
        //  If domain contains x[end]
        i0 = szz->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = szz->GetStart(2); i2 <= szz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = szz->GetStart(1); i1 <= szz->GetEnd(1); i1++)
                    sxx->SetVal(i0+ir,i1,i2) = sxx->GetVal(i0-ir+1,i1,i2);
        
        i0 = sxz->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxz->GetStart(2); i2 <= sxz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = sxz->GetStart(1); i1 <= sxz->GetEnd(1); i1++)
                    sxz->SetVal(i0+ir,i1,i2) = sxz->GetVal(i0-ir,i1,i2);
        
        i0 = sxy->GetEnd(0);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxy->GetStart(2); i2 <= sxy->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = sxy->GetStart(1); i1 <= sxy->GetEnd(1); i1++)
                    sxy->SetVal(i0+ir,i1,i2) = sxy->GetVal(i0-ir,i1,i2);
        
        //  If domain contains y[start]
        int i1 = szz->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = szz->GetStart(2); i2 <= szz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = szz->GetStart(0); i0 <= szz->GetEnd(0); i0++)
                    syy->SetVal(i0,i1-ir,i2) = syy->GetVal(i0,i1+ir,i2);
        
        i1 = syz->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = syz->GetStart(2); i2 <= syz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = syz->GetStart(0); i0 <= syz->GetEnd(0); i0++)
                    syz->SetVal(i0,i1-ir,i2) = syz->GetVal(i0,i1+ir-1,i2);
        
        i1 = sxy->GetStart(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxy->GetStart(2); i2 <= sxy->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = sxy->GetStart(0); i0 <= sxy->GetEnd(0); i0++)
                    sxy->SetVal(i0,i1-ir,i2) = sxy->GetVal(i0,i1+ir-1,i2);
        
        //  If domain contains y[end]
        i1 = szz->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = szz->GetStart(2); i2 <= szz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = szz->GetStart(0); i0 <= szz->GetEnd(0); i0++)
                    syy->SetVal(i0,i1+ir,i2) = syy->GetVal(i0,i1-ir,i2);
        
        i1 = syz->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = syz->GetStart(2); i2 <= syz->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = syz->GetStart(0); i0 <= syz->GetEnd(0); i0++)
                    syz->SetVal(i0,i1+ir,i2) = syz->GetVal(i0,i1-ir+1,i2);
        
        i1 = sxy->GetEnd(1);
#ifdef _OPENMP
#pragma omp for
#endif
        for(int i2 = sxy->GetStart(2); i2 <= sxy->GetEnd(2); i2++)
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = sxy->GetStart(0); i0 <= sxy->GetEnd(0); i0++)
                    sxy->SetVal(i0,i1+ir,i2) = sxy->GetVal(i0,i1-ir+1,i2);
        
        //  If domain contains z[start]
        if(top)
        {
            int i2 = szz->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = szz->GetStart(1)+1; i1 < szz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius; ir++)
                    for(int i0 = szz->GetStart(0); i0 <= szz->GetEnd(0); i0++)
                        szz->SetVal(i0,i1,i2-ir) = szz->GetVal(i0,i1,i2+ir);
            
            i2 = syz->GetStart(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz->GetStart(1); i1 <= syz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius + flag; ir++)
                    for(int i0 = syz->GetStart(0); i0 <= syz->GetEnd(0); i0++)
                        syz->SetVal(i0,i1,i2-ir) = (1-2*flag)*syz->GetVal(i0,i1,i2+ir-1);
            
            i2 = sxz->GetStart(2);
            //TODO: +1/-1
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = sxz->GetStart(1)+1; i1 < sxz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius + flag; ir++)
                    for(int i0 = sxz->GetStart(0)+1; i0 < sxz->GetEnd(0); i0++)
                        sxz->SetVal(i0,i1,i2-ir) =  (1-2*flag)*sxz->GetVal(i0,i1,i2+ir-1);
        }
        
        //  If domain contains z[end]
        if(bottom)
        {
            int i2 = szz->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = szz->GetStart(1)+1; i1 < szz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius; ir++)
                    for(int i0 = szz->GetStart(0); i0 <= szz->GetEnd(0); i0++)
                        szz->SetVal(i0,i1,i2+ir) = szz->GetVal(i0,i1,i2-ir);
            
            i2 = syz->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = syz->GetStart(1); i1 <= syz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius; ir++)
                    for(int i0 = syz->GetStart(0); i0 <= syz->GetEnd(0); i0++)
                        syz->SetVal(i0,i1,i2+ir) = syz->GetVal(i0,i1,i2-ir+1);
            
            i2 = sxz->GetEnd(2);
#ifdef _OPENMP
#pragma omp for
#endif
            for(int i1 = sxz->GetStart(1)+1; i1 < sxz->GetEnd(1); i1++)
                for(int ir = 1; ir < radius; ir++)
                    for(int i0 = sxz->GetStart(0)+1; i0 < sxz->GetEnd(0); i0++)
                        sxz->SetVal(i0,i1,i2+ir) = sxz->GetVal(i0,i1,i2-ir+1);
        }
    }
}
