===================================================================
ECSG - A package for anisotropic elastic wave simulations based on a stable composite staggered grid finite difference time domain numerical scheme
===================================================================

## Note: 
1. This package is able to perform both single grid and multiple (>2) grids schemes with energy-conserving condition implemented at grid refinement interface to ensure numerical stability.
2. This package relies on Madagascar/rsf interface, so install Madagascar before compiling the code.
3.

===================================================================
## Structure of the code
* **SConstruct-local/SConstruct-tacc**: python scripts for code compilation on OSX and TACC/Stampede respectively. To compile this package, copy either one of the scripts to SConstruct, and then run scons.
* **main.cpp**: driver of this program.
* **init.hpp**: set up precision macros.
* **grid.hpp**: A C struct defines grid.
* **field.hpp/field.cpp**: A C++ class defines simulation fields and related operations.
* **array.hpp/array.cpp**: A C++ class defines multi-dimensional matrices.
* **sim_utils.hpp**: list of simulation function signatures which are defined in:
  * **init_del_par_field.cpp**: functions initialize and delete parameter fields.
  * **init_del_sim.cpp**: functions initialize and delete simulation fields
  * **init_del_ghost_container.cpp**: functions initialize and delete interface ghost data.
  * **init_del_fd_coeff.cpp**: functions initialize and delete finite difference coefficients.
  * **init_src.cpp**: functions initialize a Ricker source.
  * **init_del_rec_trace**: functions initialize and delete a set of receiver traces.
* **init_del_movie.cpp**: functions initlalize and delete rsf movies.
* **init_del_res_trace.cpp**: functions initialize and 

