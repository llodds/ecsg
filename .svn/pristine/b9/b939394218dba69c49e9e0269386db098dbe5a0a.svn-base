#include "field.hpp"

void InitSim(int num_grids, int radius, bool check_energy, Grid grids[], 
             Field **&sim, Field **&vel_copy)
{
    sim = new Field* [num_grids];
    
    for(int i = 0; i < num_grids; i++)
    {
        sim[i] = new Field [9];
        
        std::string sim_strings[9] = {"vx", "vy", "vz",
            "sxx", "syy", "szz", "syz", "sxz", "sxy"};
        int sim_pos[9][3] = {
            {0,0,0},//vx
            {1,1,0},//vy
            {1,0,1},//vz
            {1,0,0},//sxx
            {1,0,0},//syy
            {1,0,0},//szz
            {1,1,1},//syz
            {0,0,1},//sxz
            {0,1,0}//sxy
        };
        
        for(int j = 0; j < 9; j++)
            sim[i][j] = Field(sim_strings[j], grids[i].grid_name,
                              sim_pos[j], grids[i].size,
                              radius);
    }
    
    if(check_energy)
    {
        vel_copy = new Field* [num_grids];
        for(int i = 0; i < num_grids; i++)
            vel_copy[i] = new Field [3];
    }
}



void DelSim(int num_grids, bool check_energy,
            Field **&sim, Field **&vel_copy)
{
    for(int i = 0; i < num_grids; i++)
        delete [] sim[i];
    delete [] sim;
    
    if(check_energy)
    {
        for(int i = 0; i < num_grids; i++)
            delete [] vel_copy[i];
        delete [] vel_copy;
    }
}