#include <rsf.hh>
#include <vector>
#include <valarray>
#include "field.hpp"
#include <iostream>

//  print grid info
void PrintGrid(Grid grid)
{
    std::cout <<"==========================="<<std::endl;
    
    std::cout << "Grid Name: ["
    <<grid.grid_name<<"]"<<std::endl;
    
    std::cout << "Grid Size(z/x/y): ["
    <<grid.size[0]<<","
    <<grid.size[1]<<","
    <<grid.size[2]<<"]"<<std::endl;
    
    std::cout << "Grid Step(z/x/y): ["
    <<grid.step[0]<<","
    <<grid.step[1]<<","
    <<grid.step[2]<<"]"<<std::endl;
    
    std::cout <<"===========================\n"<<std::endl;
}


//  Print 1D array to rsf file
void PrintToRSF(const std::vector<real> &array, std::string val_name, real dt)
{
    oRSF rsf_file(val_name.c_str());
    
    std::valarray<float> rsf_vec(array.size());
    for(int i=0; i<array.size(); i++)
        rsf_vec[i] = static_cast<float>(array[i]);
    
    rsf_file.type(SF_FLOAT);
    rsf_file.put("n1", static_cast<int>(array.size()));
    rsf_file.put("o1", static_cast<float>(0.0));
    rsf_file.put("d1", static_cast<float>(dt));
    rsf_file.put("label1", "Time");
    rsf_file.put("unit1", "ms");
    
    rsf_file.put("n2", static_cast<int>(1));
    rsf_file.put("n3", static_cast<int>(1));
    rsf_file.put("n4", static_cast<int>(1));
    
    rsf_file << rsf_vec;
}


//  Print Field movie to rsf file
void PrintToRSF(const Field *movie, int nframe,
                Grid grid, real time_interval,
                std::string movie_name)
{
    oRSF rsf_file(movie_name.c_str());
    
    rsf_file.type(SF_FLOAT);
    
    rsf_file.put("n1", static_cast<int>(movie[0].GetEnd(0)-movie[0].GetStart(0)+1));
    rsf_file.put("n2", static_cast<int>(movie[0].GetEnd(1)-movie[0].GetStart(1)+1));
    rsf_file.put("n3", static_cast<int>(movie[0].GetEnd(2)-movie[0].GetStart(2)+1));
    rsf_file.put("n4", static_cast<int>(nframe));
    
    rsf_file.put("d1", static_cast<float>(grid.step[0]));
    rsf_file.put("d2", static_cast<float>(grid.step[1]));
    rsf_file.put("d3", static_cast<float>(grid.step[2]));
    rsf_file.put("d4", static_cast<float>(time_interval));
    
    rsf_file.put("o1", static_cast<float>(grid.origin[0]));
    rsf_file.put("o2", static_cast<float>(grid.origin[1]));
    rsf_file.put("o3", static_cast<float>(grid.origin[2]));
    rsf_file.put("o4", static_cast<float>(0)); //assume the movie starts from 0 ms.
    
    rsf_file.put("label1", "Distance_Z");rsf_file.put("unit1", "m");
    rsf_file.put("label2", "Distance_X");rsf_file.put("unit2", "m");
    rsf_file.put("label3", "Distance_Y");rsf_file.put("unit3", "m");
    rsf_file.put("label4", "Time");rsf_file.put("unit4", "ms");
    
    std::valarray<float> trace(movie[0].GetEnd(0)-movie[0].GetStart(0)+1);
    for(int it = 0; it < nframe; it++)
        for(int i2 = movie[0].GetStart(2); i2 <= movie[0].GetEnd(2); i2++)
            for(int i1 = movie[0].GetStart(1); i1 <= movie[0].GetEnd(1); i1++){
                for(int i0 = movie[0].GetStart(0); i0 <= movie[0].GetEnd(0); i0++)
                    trace[i0-movie[0].GetStart(0)]=static_cast<float>(movie[it].GetVal(i0,i1,i2));
                rsf_file << trace;
            }
}


//print a single field to existing oRSF
void PrintFieldToRSF(const Field &f, oRSF *out)
{
    
    int size0 = f.GetEnd(0)-f.GetStart(0)+1;
    
    std::valarray<float> trace(size0);
    
    for(int i2 = f.GetStart(2); i2 <= f.GetEnd(2); i2++)
        for(int i1 = f.GetStart(1); i1 <= f.GetEnd(1); i1++)
        {
            for(int i0 = f.GetStart(0); i0 <= f.GetEnd(0); i0++)
                trace[i0-f.GetStart(0)] = f.GetVal(i0,i1,i2);
            
            *out << trace;
        }
}


//print receiver data samples
void PrintRecData(//input
                  int n_rec,
                  std::vector<real>* rec_data,
                  real dt,
                  //output
                  oRSF &rec_out)
{
    rec_out.type(SF_FLOAT);
    rec_out.put("n1", static_cast<int>(rec_data[0].size()));
    rec_out.put("o1", static_cast<float>(0));
    rec_out.put("d1", static_cast<float>(dt));
    rec_out.put("label1", "Time"); rec_out.put("unit1","ms");
    
    rec_out.put("n2", static_cast<int>(n_rec));
    
    std::valarray<float> trace(rec_data[0].size());
    for(int j = 0; j < n_rec; j++)
    {
        for(int i = 0; i < rec_data[0].size(); i++)
            trace[i] = rec_data[j][i];
        rec_out << trace;
    }
}
