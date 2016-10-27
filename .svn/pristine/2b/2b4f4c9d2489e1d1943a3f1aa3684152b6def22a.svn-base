#include "field.hpp"
#include <cassert>
#include "sim_utils.hpp"

//setup material parameter field
void InitParField(//input
                  int num_grids,
                  bool check_energy,
                  //output
                  Grid *&grids,
                  Field *&buoys,
                  Field **&stiffs_obj,
                  Field ***&stiffs_ptr,
                  Field **&densities_obj,
                  Field ***&densities_ptr,
                  Field **&compliances_obj,
                  Field ***&compliances_ptr)
{
    //=====================================================
    //init grids/buoys by reading in buoy files
    grids = new Grid [num_grids];
    buoys = new Field [num_grids];
    for(int i = 1; i <= num_grids; i++)
    {
        std::string buoy_name = "buoy_" + std::to_string(static_cast<long long>(i));
        iRSF in_buoy(buoy_name.c_str());
        
        int size[3];
        real step[3];
        
        for(int j = 0; j < 3; j++)
        {
            std::string size_name = "n" + std::to_string(static_cast<long long>(j+1));
            in_buoy.get(size_name.c_str(), size[j], -1);
            size[j] = size[j] - 1;
            assert(size[j] > 0);
            
            std::string step_name = "d" + std::to_string(static_cast<long long>(j+1));
            in_buoy.get(step_name.c_str(), step[j], -1);
            assert(step[j] > 0);
        }
        
        
        //init grids
        //mesh size: grids[0] < grids[1] < grid[2] < ...
        Grid& current_grid = grids[i-1];
        current_grid.grid_name = "grid_" + std::to_string(static_cast<long long>(i));
        //so grids names are "grid1", "grid2", "grid3", ...
        for(int j = 0; j < 3; j++)
        {
            current_grid.size[j] = size[j];
            current_grid.step[j] = step[j];
            current_grid.origin[j] = 0;
        }
        if(i >= 2)
            current_grid.origin[2] = grids[i-2].origin[2]
            + grids[i-2].size[2]*grids[i-2].step[2];
        
        //init buoy fields
        int buoy_pos[3] = {0,0,0}; //buoy is always placed on the primary grid
        buoys[i-1] = Field(buoy_name, current_grid.grid_name,
                           buoy_pos, current_grid.size,
                           0);
        std::valarray<float> trace(current_grid.size[0]+1);
        for(int iz=0; iz<=current_grid.size[2]; iz++)
            for(int iy=0; iy<=current_grid.size[1]; iy++)
            {
                in_buoy >> trace;
                for(int ix=0; ix<=current_grid.size[0]; ix++)
                    buoys[i-1].SetVal(ix,iy,iz) = trace[ix];
            }
    }
    
    
    //=====================================
    //init stiffs by reading in stiff files
    //TODO: add other anisotropic mode
    stiffs_obj = new Field* [num_grids]; //store the Field obj
    stiffs_ptr = new Field** [num_grids]; //store the Field ptr, because
    //multiple Field ptrs could point to the same Field obj
    std::string stiff_name[5] = {"s11_", "s13_", "s44_", "s55_", "s66_"};
    int stiff_pos[5][3] = {
        {1,0,0},//s11
        {1,0,0},//s13
        {1,1,1},//s44
        {0,0,1},//s55
        {0,1,0}//s66
    };
    
    for(int i = 1; i <= num_grids; i++)
    {
        stiffs_obj[i-1] = new Field [5];
        //stiffs_obj[][0]: s11
        //stiffs_obj[][1]: s13
        //stiffs_obj[][2]: s44
        //stiffs_obj[][3]: s55
        //stiffs_obj[][4]: s66
        
        for(int j = 0; j < 2; j++)
        {
            std::string s_name = stiff_name[j] + std::to_string(static_cast<long long>(i));
            iRSF in_s(s_name.c_str());
            
            int primary_pos[3] = {0,0,0};
            Grid& current_grid = grids[i-1];
            Field tmp_s(s_name, current_grid.grid_name,
                        primary_pos, current_grid.size,
                        0);
            std::valarray<float> trace(current_grid.size[0]+1);
            for(int iz=0; iz<=current_grid.size[2]; iz++)
                for(int iy=0; iy<=current_grid.size[1]; iy++)
                {
                    in_s >> trace;
                    for(int ix=0; ix<=current_grid.size[0]; ix++)
                        tmp_s.SetVal(ix,iy,iz) = trace[ix];
                }
            stiffs_obj[i-1][j] = tmp_s;
        }
        stiffs_obj[i-1][2] = stiffs_obj[i-1][0]-stiffs_obj[i-1][1];
        
        for(int j = 0; j < 2; j++)
            stiffs_obj[i-1][j] = ReshapeField(stiffs_obj[i-1][j], stiff_pos[j]);
        for(int j = 3; j < 5; j++)
            stiffs_obj[i-1][j] = ReshapeField(stiffs_obj[i-1][2], stiff_pos[j], stiff_name[j] + std::to_string(static_cast<long long>(i)));
        stiffs_obj[i-1][2] = ReshapeField(stiffs_obj[i-1][2], stiff_pos[2], stiff_name[2] + std::to_string(static_cast<long long>(i)));
        
        //assign ptrs
        stiffs_ptr[i-1] = new Field* [9];
        //stiffs_ptr[][0]: s11
        //stiffs_ptr[][0]: s22
        //stiffs_ptr[][0]: s33
        //stiffs_ptr[][0]: s23
        //stiffs_ptr[][0]: s13
        //stiffs_ptr[][0]: s12
        //stiffs_ptr[][0]: s44
        //stiffs_ptr[][0]: s55
        //stiffs_ptr[][0]: s66
        
        for(int j = 0; j < 3; j++)
        {
            stiffs_ptr[i-1][j] = stiffs_obj[i-1];
        }
        for(int j = 3; j < 6; j++)
        {
            stiffs_ptr[i-1][j] = stiffs_obj[i-1]+1;
        }
        stiffs_ptr[i-1][6] = stiffs_obj[i-1]+2;
        stiffs_ptr[i-1][7] = stiffs_obj[i-1]+3;
        stiffs_ptr[i-1][8] = stiffs_obj[i-1]+4;
    }
    
    
    //=====================================================
    //init other parameter fields for computing the energy
    if(check_energy)
    {
        densities_ptr = new Field** [num_grids];
        densities_obj = new Field* [num_grids];
        
        compliances_ptr = new Field** [num_grids];
        compliances_obj = new Field* [num_grids];
        
        for(int i = 0; i < num_grids; i++)
        {
            densities_ptr[i] = new Field* [3];
            densities_obj[i] = new Field [3];
            
            compliances_ptr[i] = new Field* [9];
            compliances_obj[i] = new Field [9];
        }
        
        for(int i = 0; i < num_grids; i++)
            InvertParameters(densities_ptr[i], densities_obj[i],
                             compliances_ptr[i], compliances_obj[i],
                             buoys[i], stiffs_ptr[i],
                             ISOTROPIC_MODE);
    }
}

//delete material parameter field
void DelParField(//input
                    int num_grids,
                    bool check_energy,
                    //output
                    Grid *&grids,
                    Field *&buoys,
                    Field **&stiffs_obj,
                    Field ***&stiffs_ptr,
                    Field **&densities_obj,
                    Field ***&densities_ptr,
                    Field **&compliances_obj,
                    Field ***&compliances_ptr)
{
    delete [] grids;
    
    delete [] buoys;
    
    for(int i = 0; i < num_grids; i++)
    {
        delete [] stiffs_obj[i];
        delete [] stiffs_ptr[i];
    }
    delete [] stiffs_obj;
    delete [] stiffs_ptr;
    
    if(check_energy)
    {
        for(int i = 0; i < num_grids; i++)
        {
            delete [] densities_obj[i];
            delete [] densities_ptr[i];
            delete [] compliances_obj[i];
            delete [] compliances_ptr[i];
        }
        delete [] densities_ptr;
        delete [] densities_obj;
        delete [] compliances_ptr;
        delete [] compliances_obj;
    }
}



void InvertParameters(Field *density_ptr[3], Field density_obj[3],
                      Field *compliance_ptr[9], Field compliance_obj[9],
                      const Field &buoy, Field *stiffness_ptr[9], int aniso_mode)
{
    
    //=== invert for density ===
    //1. compute density for vx
    density_obj[0] = InvField(buoy);
    density_obj[0].SetFieldName("buoy_x");
    density_ptr[0] = &density_obj[0];
    
    //2. compute density for vy
    int pos_vy[3] = {1,1,0};
    density_obj[1] = ReshapeField(buoy,pos_vy);
    density_obj[1] = InvField(density_obj[1]);
    density_obj[1].SetFieldName("buoy_y");
    density_ptr[1] = &density_obj[1];
    
    //3. compute density for vz
    int pos_vz[3] = {1,0,1};
    density_obj[2] = ReshapeField(buoy,pos_vz);
    density_obj[2] = InvField(density_obj[2]);
    density_obj[2].SetFieldName("buoy_z");
    density_ptr[2] = &density_obj[2];
    
    
    if(aniso_mode == ISOTROPIC_MODE)
    {
        //=== invert for compliance ===
        Field* s11 = stiffness_ptr[0];
        Field* s13 = stiffness_ptr[4];
        Field* s44 = stiffness_ptr[6];
        Field* s55 = stiffness_ptr[7];
        Field* s66 = stiffness_ptr[8];
        
        Field base = (*s11)*(*s13)*(*s13)*3
        - (*s13)*(*s13)*(*s13)*2
        - (*s11)*(*s11)*(*s11);
        base = InvField(base);
        
        //1. compute for c11
        compliance_obj[0] = ((*s13)*(*s13) - (*s11)*(*s11))*base;
        compliance_obj[0].SetFieldName("c11");
        //2. compute for c13
        compliance_obj[4] = ((*s13)*(*s11) - (*s13)*(*s13))*base;
        compliance_obj[4].SetFieldName("c13");
        //3. compute for c44
        compliance_obj[6] = InvField(*s44);
        compliance_obj[6].SetFieldName("c44");
        //4. compute for c55
        compliance_obj[7] = InvField(*s55);
        compliance_obj[7].SetFieldName("c55");
        //5. compute for c66
        compliance_obj[8] = InvField(*s66);
        compliance_obj[8].SetFieldName("c66");
        
        
        for(int i=0; i<3; i++) compliance_ptr[i] = &compliance_obj[0];
        for(int i=3; i<6; i++) compliance_ptr[i] = &compliance_obj[4];
        for(int i=6; i<9; i++) compliance_ptr[i] = &compliance_obj[i];
    }
    else
    {
        std::cerr << "ERROR: Non-isotropic cases haven't been implemented yet" << std::endl;
        exit(-1);
    }
}
