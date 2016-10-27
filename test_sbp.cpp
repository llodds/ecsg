#include "field.hpp"
#include <iostream>
#include "init.hpp"

void TestSBP(Field** sim)
{
    //test1: (u, D[-,z]sxz)_u + (D[+,z]u, sxz)_sxz = 0
    Field *vx = &(sim[0][0]);
    Field *sxz = &(sim[0][7]);
    real SBP_res1 = 0;
    for(int iz = vx->GetStart(2); iz <= vx->GetEnd(2); iz++)
        for(int iy = vx->GetStart(1); iy <= vx->GetEnd(1); iy++)
            for(int ix = vx->GetStart(0); ix <= vx->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==vx->GetStart(2)))
                *(1-0.5*(iz==vx->GetEnd(2)))
                *(1-0.5*(iy==vx->GetStart(1)))
                *(1-0.5*(iy==vx->GetEnd(1)))
                *(1-0.5*(ix==vx->GetStart(0)))
                *(1-0.5*(ix==vx->GetEnd(0)));
                //                real mul = (iz!=vx->GetStart(2))
                //                *(iz!=vx->GetEnd(2))
                //                *(iy!=vx->GetStart(1))
                //                *(iy!=vx->GetEnd(1))
                //                *(ix!=vx->GetStart(0))
                //                *(ix!=vx->GetEnd(0));
                SBP_res1 += mul*vx->GetVal(ix,iy,iz)*(sxz->GetVal(ix,iy,iz)-sxz->GetVal(ix,iy,iz-1));
            }
    for(int iz = sxz->GetStart(2); iz <= sxz->GetEnd(2); iz++)
        for(int iy = sxz->GetStart(1); iy <= sxz->GetEnd(1); iy++)
            for(int ix = sxz->GetStart(0); ix <= sxz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==sxz->GetStart(1)))
                *(1-0.5*(iy==sxz->GetEnd(1)))
                *(1-0.5*(ix==sxz->GetStart(0)))
                *(1-0.5*(ix==sxz->GetEnd(0)));
                SBP_res1 += mul*(vx->GetVal(ix,iy,iz+1)-vx->GetVal(ix,iy,iz))*sxz->GetVal(ix,iy,iz);
            }
    std::cerr << "test1: (u, D[-,z]sxz)_u + (D[+,z]u, sxz)_sxz = 0" <<std::endl;
    std::cerr << KRED << "\t\t SBP_res1 = " << SBP_res1 << RST << std::endl;
    
    
    //test2: (u, D[-,x]sxx)_u + (D[+,x]u, sxx)_sxx = 0
    //Field *vx = &(sim[0][0]);
    Field *sxx = &(sim[0][3]);
    real SBP_res2 = 0;
    for(int iz = vx->GetStart(2); iz <= vx->GetEnd(2); iz++)
        for(int iy = vx->GetStart(1); iy <= vx->GetEnd(1); iy++)
            for(int ix = vx->GetStart(0); ix <= vx->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==vx->GetStart(2)))
                *(1-0.5*(iz==vx->GetEnd(2)))
                *(1-0.5*(iy==vx->GetStart(1)))
                *(1-0.5*(iy==vx->GetEnd(1)))
                *(1-0.5*(ix==vx->GetStart(0)))
                *(1-0.5*(ix==vx->GetEnd(0)));
                SBP_res2 += mul*vx->GetVal(ix,iy,iz)*(sxx->GetVal(ix,iy,iz)-sxx->GetVal(ix-1,iy,iz));
            }
    for(int iz = sxx->GetStart(2); iz <= sxx->GetEnd(2); iz++)
        for(int iy = sxx->GetStart(1); iy <= sxx->GetEnd(1); iy++)
            for(int ix = sxx->GetStart(0); ix <= sxx->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==sxx->GetStart(2)))
                *(1-0.5*(iz==sxx->GetEnd(2)))
                *(1-0.5*(iy==sxx->GetStart(1)))
                *(1-0.5*(iy==sxx->GetEnd(1)));
                SBP_res2 += mul*(vx->GetVal(ix+1,iy,iz)-vx->GetVal(ix,iy,iz))*sxx->GetVal(ix,iy,iz);
            }
    std::cerr << "test2: (u, D[-,x]sxx)_u + (D[+,x]u, sxx)_sxx = 0" <<std::endl;
    std::cerr << KRED << "\t\t SBP_res2 = " << SBP_res2 << RST << std::endl;
    
    
    
    //test3: (u, D[-,x]sxy)_u + (D[+,x]u, sxy)_sxy = 0
    //Field *vx = &(sim[0][0]);
    Field *sxy = &(sim[0][8]);
    real SBP_res3 = 0;
    for(int iz = vx->GetStart(2); iz <= vx->GetEnd(2); iz++)
        for(int iy = vx->GetStart(1); iy <= vx->GetEnd(1); iy++)
            for(int ix = vx->GetStart(0); ix <= vx->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==vx->GetStart(2)))
                *(1-0.5*(iz==vx->GetEnd(2)))
                *(1-0.5*(iy==vx->GetStart(1)))
                *(1-0.5*(iy==vx->GetEnd(1)))
                *(1-0.5*(ix==vx->GetStart(0)))
                *(1-0.5*(ix==vx->GetEnd(0)));
                SBP_res3 += mul*vx->GetVal(ix,iy,iz)*(sxy->GetVal(ix,iy,iz)-sxy->GetVal(ix,iy-1,iz));
            }
    for(int iz = sxy->GetStart(2); iz <= sxy->GetEnd(2); iz++)
        for(int iy = sxy->GetStart(1); iy <= sxy->GetEnd(1); iy++)
            for(int ix = sxy->GetStart(0); ix <= sxy->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==sxy->GetStart(2)))
                *(1-0.5*(iz==sxy->GetEnd(2)))
                *(1-0.5*(ix==sxy->GetStart(0)))
                *(1-0.5*(ix==sxy->GetEnd(0)));
                SBP_res3 += mul*(vx->GetVal(ix,iy+1,iz)-vx->GetVal(ix,iy,iz))*sxy->GetVal(ix,iy,iz);
            }
    std::cerr << "test3: (u, D[-,x]sxy)_u + (D[+,x]u, sxy)_sxy = 0" <<std::endl;
    std::cerr <<  KRED << "\t\t SBP_res3 = " << SBP_res3 << RST << std::endl;
    
    
    
    
    
    //test4: (v, D[-,z]syz)_v + (D[+,z]v, syz)_syz = 0
    Field *vy = &(sim[0][1]);
    Field *syz = &(sim[0][6]);
    real SBP_res4 = 0;
    for(int iz = vy->GetStart(2); iz <= vy->GetEnd(2); iz++)
        for(int iy = vy->GetStart(1); iy <= vy->GetEnd(1); iy++)
            for(int ix = vy->GetStart(0); ix <= vy->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iz==vy->GetStart(2)))
                *(1-0.5*(iz==vy->GetEnd(2)));
                SBP_res4 += mul*vy->GetVal(ix,iy,iz)*(syz->GetVal(ix,iy,iz)-syz->GetVal(ix,iy,iz-1));
            }
    for(int iz = syz->GetStart(2); iz <= syz->GetEnd(2); iz++)
        for(int iy = syz->GetStart(1); iy <= syz->GetEnd(1); iy++)
            for(int ix = syz->GetStart(0); ix <= syz->GetEnd(0); ix++)
            {
                SBP_res4 += (vy->GetVal(ix,iy,iz+1)-vy->GetVal(ix,iy,iz))*syz->GetVal(ix,iy,iz);
            }
    std::cerr << "test4: (v, D[-,z]syz)_v + (D[+,z]v, syz)_syz = 0" <<std::endl;
    std::cerr <<  KRED << "\t\t SBP_res4 = " << SBP_res4 << RST << std::endl;
    
    
    
    
    //test5: (w, D[+,z]szz)_w + (D[-,z]w, szz)_szz = 0
    Field *vz = &(sim[0][2]);
    Field *szz = &(sim[0][5]);
    real SBP_res5 = 0;
    for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
        for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
            for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==vz->GetStart(1)))
                *(1-0.5*(iy==vz->GetEnd(1)));
                SBP_res5 += mul*vz->GetVal(ix,iy,iz)*(szz->GetVal(ix,iy,iz+1)-szz->GetVal(ix,iy,iz));
            }
    for(int iz = szz->GetStart(2); iz <= szz->GetEnd(2); iz++)
        for(int iy = szz->GetStart(1); iy <= szz->GetEnd(1); iy++)
            for(int ix = szz->GetStart(0); ix <= szz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==szz->GetStart(1)))
                *(1-0.5*(iy==szz->GetEnd(1)))
                *(1-0.5*(iz==szz->GetStart(2)))
                *(1-0.5*(iz==szz->GetEnd(2)));
                SBP_res5 += mul*(vz->GetVal(ix,iy,iz)-vz->GetVal(ix,iy,iz-1))*szz->GetVal(ix,iy,iz);
            }
    std::cerr << "test5: (w, D[+,z]szz)_w + (D[-,z]w, szz)_szz = 0" <<std::endl;
    std::cerr <<  KRED << "\t\t SBP_res5 = " << SBP_res5 << RST << std::endl;
    
    
    
    //test6: (w, D[+,x]sxz)_w + (D[-,x]w, sxz)_sxz = 0
    real SBP_res6 = 0;
    for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
        for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
            for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==vz->GetStart(1)))
                *(1-0.5*(iy==vz->GetEnd(1)));
                SBP_res6 += mul*vz->GetVal(ix,iy,iz)*(sxz->GetVal(ix+1,iy,iz)-sxz->GetVal(ix,iy,iz));
            }
    for(int iz = sxz->GetStart(2); iz <= sxz->GetEnd(2); iz++)
        for(int iy = sxz->GetStart(1); iy <= sxz->GetEnd(1); iy++)
            for(int ix = sxz->GetStart(0); ix <= sxz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==sxz->GetStart(1)))
                *(1-0.5*(iy==sxz->GetEnd(1)))
                *(1-0.5*(ix==sxz->GetStart(0)))
                *(1-0.5*(ix==sxz->GetEnd(0)));
                SBP_res6 += mul*(vz->GetVal(ix,iy,iz)-vz->GetVal(ix-1,iy,iz))*sxz->GetVal(ix,iy,iz);
            }
    std::cerr << "test6: (w, D[+,x]sxz)_w + (D[-,x]w, sxz)_sxz = 0" <<std::endl;
    std::cerr <<  KRED << "\t\t SBP_res6 = " << SBP_res6 << RST << std::endl;
    
    
    //test7: (w, D[+,y]syz)_w + (D[-,y]w, syz)_syz = 0
    real SBP_res7 = 0;
    for(int iz = vz->GetStart(2); iz <= vz->GetEnd(2); iz++)
        for(int iy = vz->GetStart(1); iy <= vz->GetEnd(1); iy++)
            for(int ix = vz->GetStart(0); ix <= vz->GetEnd(0); ix++)
            {
                real mul = (1-0.5*(iy==vz->GetStart(1)))
                *(1-0.5*(iy==vz->GetEnd(1)));
                SBP_res7 += mul*vz->GetVal(ix,iy,iz)*(syz->GetVal(ix,iy+1,iz)-syz->GetVal(ix,iy,iz));
            }
    for(int iz = syz->GetStart(2); iz <= syz->GetEnd(2); iz++)
        for(int iy = syz->GetStart(1); iy <= syz->GetEnd(1); iy++)
            for(int ix = syz->GetStart(0); ix <= syz->GetEnd(0); ix++)
            {
                SBP_res7 += (vz->GetVal(ix,iy,iz)-vz->GetVal(ix,iy-1,iz))*syz->GetVal(ix,iy,iz);
            }
    std::cerr << "test6: (w, D[+,y]syz)_w + (D[-,y]w, syz)_syz = 0" <<std::endl;
    std::cerr <<  KRED << "\t\t SBP_res7 = " << SBP_res7 << RST << std::endl;
}