#include "init.hpp"
#include "field.hpp"

/*
 sim_ptr[0]: vx
 sim_ptr[1]: vy
 sim_ptr[2]: vz
 sim_ptr[3]: sxx
 sim_ptr[4]: syy
 sim_ptr[5]: szz
 sim_ptr[6]: syz
 sim_ptr[7]: sxz
 sim_ptr[8]: sxy
 
 vel_copy[0]: vx_copy
 vel_copy[1]: vy_copy
 vel_copy[2]: vz_copy
 
 density_ptr[0]: density_vx
 density_ptr[1]: density_vy
 density_ptr[2]: density_vz
 
 compliance_ptr[0]: c11
 compliance_ptr[1]: c22
 compliance_ptr[2]: c33
 compliance_ptr[3]: c23
 compliance_ptr[4]: c13
 compliance_ptr[5]: c12
 compliance_ptr[6]: c44
 compliance_ptr[7]: c55
 compliance_ptr[8]: c66
 */

real ComputeEnergy(const Field *sim, const Field vel_copy[3], Field * density_ptr[3], Field *compliance_ptr[9], const real grid_step[3])
{
    real res = 0;
    //add density*(vx|vy|vz)^2
    for(int i = 0; i < 3; i++)
        res += TriMulFields(*(density_ptr[i]), sim[i], vel_copy[i]);
    //add compliance*(sxx|syy|szz)^2
    for(int i = 0; i < 3; i++)
        res += TriMulFields(*(compliance_ptr[i]), sim[i+3], sim[i+3]);
    //add compliance*(syy*szz|sxx*szz|sxx*syy)^2
    res += 2*TriMulFields(*(compliance_ptr[3]), sim[4], sim[5]);
    res += 2*TriMulFields(*(compliance_ptr[4]), sim[3], sim[5]);
    res += 2*TriMulFields(*(compliance_ptr[5]), sim[3], sim[4]);
    //add compliance*(syz|sxz|sxy)^2
    for(int i = 6; i < 9; i++)
        res += 2*TriMulFields(*(compliance_ptr[i]), sim[i], sim[i]);
    
    //finally, multiply by the half cell size
    res *= 0.5*grid_step[2]*grid_step[1]*grid_step[0];
    return res;
}
