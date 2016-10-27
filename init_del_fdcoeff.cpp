#include "init.hpp"
#include <vector>
#include <iostream>
#include "grid.hpp"
#include "sim_utils.hpp"

void InitFDcoeff(int num_grids, Grid grids[], int radius, real dt, real ***& coeff)
{
    coeff = new real** [num_grids];
    
    std::vector<real> fd_coeff=FDCoeff(radius);
    for(int i = 0; i < num_grids; i++)
    {
        coeff[i] = new real* [3];
        for(int j = 0; j < 3; j++)
        {
            coeff[i][j] = new real [radius];
            for(int k = 0; k < radius; k++)
                coeff[i][j][k] = fd_coeff[k]*dt/grids[i].step[j];
        }
    }
}


void DelFDCoeff(int num_grids, real ***& coeff)
{
    for(int i = 0; i < num_grids; i++)
    {
        for(int j = 0; j < 3; j++)
            delete [] coeff[i][j];
        delete [] coeff[i];
    }
    delete [] coeff;
}


//TODO: test if the precison changes if using MKL to solve the LES
std::vector<real> FDCoeff(int radius)
{
    std::vector<real> res(radius, 0);
    
    switch (radius) {
        case 1:
            res[0] = (real)1;
            break;
        case 2:
            res[0] = (real)9/8;
            res[1] = -(real)1/24;
            break;
        case 3:
            res[0] = (real)75/64;
            res[1] = -(real)25/384;
            res[2] = (real)3/640;
            break;
        case 4:
            res[0] = (real)1225/1024;
            res[1] = -(real)245/3072;
            res[2] = (real)49/5120;
            res[3] = -(real)5/7168;
            break;
        case 5:
            res[0] = (real)711/587;
            res[1] = -(real)158/1761;
            res[2] = (real)296/21383;
            res[3] = -(real)25/14159;
            res[4] = (real)17/143243;
            break;
        case 6:
            res[0] = (real)1225/1003;
            res[1] = -(real)338/3487;
            res[2] = (real)35/2006;
            res[3] = -(real)127/42800;
            res[4] = (real)19/52924;
            res[5] = -(real)6/274627;
            break;
        default:
        {
            std::cerr<< "ERROR: FD coefficents of radius="<<radius
            <<" have not been computed."<<std::endl;
            std::cerr << "File [" << __FILE__ << "], Line [" << __LINE__ << "]." << std::endl;
            exit(-1);
        }
    }
    
    return res;
}