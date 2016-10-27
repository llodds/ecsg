===================================================================
ECSG - A package for anisotropic elastic wave simulations based on a stable composite staggered grid finite difference time domain numerical scheme
===================================================================

## Note: 
##1. This code is able to perform both single grid and multiple (>2) grids schemes with energy-conserving condition implemented at grid refinement interface to ensure numerical stability.
##2. This code is able to 

===================================================================
## Structure of the code
* **SConstruct-local/SConstruct-tacc**: python scripts for code compilation on OSX and TACC/Stampede respectively. To compile this package, copy either one of the scripts to SConstruct, and then run scons.
* **main.cpp**: driver of this program.
* **init.hpp**: set up precision.
* **grid.hpp**: A C struct defines grid.
* **field.hpp/field.cpp**: A C++ class defines simulation fields and related operations.
## array.hpp/array.cpp: A C++ class defines multi-dimensional matrices.
## sim_utils.hpp: 

init.hpp: contains pre-defined macros

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
