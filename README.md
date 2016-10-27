===================================================================
ECSG -- A package for Anisotropic Elastic Wave Simulations Based on A Stable Composite Staggered Grid Finite Difference Time Domain Numerical Scheme
===================================================================

## Note: 
1. this code is able to perform both single grid and multiple grids (with energy-conserving interface
2. 

===================================================================

  * init.hpp: contains pre-defined macros

  * grid.hpp: contains the grid structure
  
  * field.hpp : contains the field class and field-related non-member functions (field.cpp)
    
  * sim_utils.hpp : contains a list of functions to assist simulations
     (add_src.cpp,
       compute_energy.cpp,
       invert_parameters.cpp,
       update_stress.cpp,
       update_vel.cpp,
       print_grid.cpp,
       ricker.cpp,
       fd_coeff.cpp)
    
  * main.cpp: main program
    
  * esg/test/SConstruct: contains a 3D test.
