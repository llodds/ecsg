#include "field.hpp"
#include <rsf.hh>
#include "sim_utils.hpp"

void InitMovie(int num_grids,
               bool movie_flags[3],
               Field **sim,
               Grid grids[],
               int movie_nt,
               real movie_dt,
               oRSF*** &movie)
{
    std::string vel_name[3] = {"vx", "vy", "vz"};
    
    for(int j=0; j<3; j++)
        if(movie_flags[j])
        {
            movie[j] = new oRSF* [num_grids];
            for(int i=0; i<num_grids; i++)
            {
                std::string movie_name = "movie_" + vel_name[j] + "_" + std::to_string(static_cast<long long>(i+1));
                movie[j][i] = new oRSF(movie_name.c_str());
                InitRSFMovie(sim[i][j],movie_nt,grids[i],movie_dt,movie[j][i]);
            }
        }
}


//init oRSF
void InitRSFMovie(const Field &f, int nframe, Grid grid, real time_interval, oRSF *rsf_file)
{
    rsf_file->type(SF_FLOAT);
    
    rsf_file->put("n1", static_cast<int>(f.GetEnd(0)-f.GetStart(0)+1));
    rsf_file->put("n2", static_cast<int>(f.GetEnd(1)-f.GetStart(1)+1));
    rsf_file->put("n3", static_cast<int>(f.GetEnd(2)-f.GetStart(2)+1));
    rsf_file->put("n4", static_cast<int>(nframe));
    
    rsf_file->put("d1", static_cast<float>(grid.step[0]));
    rsf_file->put("d2", static_cast<float>(grid.step[1]));
    rsf_file->put("d3", static_cast<float>(grid.step[2]));
    rsf_file->put("d4", static_cast<float>(time_interval));
    
    rsf_file->put("o1", static_cast<float>(grid.origin[0]));
    rsf_file->put("o2", static_cast<float>(grid.origin[1]));
    rsf_file->put("o3", static_cast<float>(grid.origin[2]));
    rsf_file->put("o4", static_cast<float>(0)); //assume the movie starts from 0 ms.
    
    rsf_file->put("label1", "X");rsf_file->put("unit1", "m");
    rsf_file->put("label2", "Y");rsf_file->put("unit2", "m");
    rsf_file->put("label3", "Z");rsf_file->put("unit3", "m");
    rsf_file->put("label4", "Time");rsf_file->put("unit4", "ms");
}


void DelMovie(int num_grids,
              bool movie_flags[3],
              oRSF ***movie)
{
    for(int j = 0; j < 3; j++)
        if(movie_flags[j])
        {
            for(int i=0; i<num_grids; i++)
                delete movie[j][i];
            delete [] movie[j];
        }
    delete [] movie;
}

void PrintMovie(int num_grids,
                bool movie_flags[3],
                Field **sim,
                oRSF ***movie)
{
    for(int j = 0; j < 3; j++)
        if(movie_flags[j])
        {
            for(int i=0; i<num_grids; i++)
                PrintFieldToRSF(sim[i][j], movie[j][i]);
        }
}