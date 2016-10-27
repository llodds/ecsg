#include <vector>

void InitRes(int num_grids,
             std::vector<float> **&res_vec)
{
    if(num_grids>1)
        res_vec = new std::vector<float>* [num_grids-1];
    for(int i = 0; i < num_grids-1; i++)
        res_vec[i] = new std::vector<float> [9];
}


void DelRes(int num_grids,
            std::vector<float> **&res_vec)
{
    for(int i = 0; i < num_grids-1; i++)
        delete [] res_vec[i];
    delete [] res_vec;
}