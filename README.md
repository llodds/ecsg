===================================================================
#ECSG - A package for 3D anisotropic elastic wave simulations based on a stable composite staggered grid finite difference time domain numerical scheme

===================================================================

### Why composite grid?
It is well known that finite difference scheme suffers from grid dispersion. To minimize this numerical aritifact, we need to sample each wavelength by at least a certain number of grid points. When the simulation domain contains a low velocity zone, for example the shallow region, using the classical uniform grid scheme (same mesh size over the entire domain) would require a very fine mesh. In contrast, a composite grid scheme uses fine mesh for low velocity region and coarse mesh for high velocity region, so it needs less grid points (less memory) to simulate the same region. 

### Why we need stability?
In numerical analysis, numerical stability is always an important property we are looking for when we design a numerical scheme. Stability means a perturbation in the initial solution (e.g., machine error) will not blow up over time. Lax-Richtmyer theorem claims that for a linear consistent numerical scheme (e.g., finite difference scheme based on the wave equation), stability is a sufficient and necessary condition for convergence, which measures how close the numerical solution will be to the true solution if we keep reducing the mesh size. Because of these reasons, we need a stable numerical scheme.

### What is energy method and how we use the energy method to achieve numerical stability for this composite staggered grid FDTD scheme based on the elastic wave equation?
After figuring out how to formulate the energy for the numerical solution on a uniform grid, we compute the energy on the composite grid by summing up energies on each individual uniform grid, 


===================================================================

## Note: 
1. This package is able to perform both single grid and multiple (>=2) grids schemes with energy-conserving condition implemented at grid refinement interfaces to ensure numerical stability.
2. This package relies on [Madagascar](https://github.com/ahay/src), OpenMP and [scons](http://scons.org), so have them installed or loaded before compiling the code.
3. Currently the scheme is based on FDTD of second order in time and space. Free surface boundary condition is implemented on top, and rigid (Dirichlet) boundary condition is implemented on other three sides of the rectangular domain.

===================================================================

## Structure
* **SConstruct-local/SConstruct-tacc**: python scripts for code compilation on OSX and TACC/Stampede respectively. To compile this package, copy either one of the scripts to SConstruct, and then run scons.
* **main.cpp**: driver of this program.
* **init.hpp**: set up precision macros.
* **grid.hpp**: A C struct defines grid.
* **field.hpp/field.cpp**: A C++ class defines simulation fields and related operations.
* **array.hpp/array.cpp**: A C++ class defines multi-dimensional matrices.
* **sim_utils.hpp**: list of simulation function signatures with functions defined in:
  * **init_del_par_field.cpp**: functions initialize and delete parameter fields.
  * **init_del_sim.cpp**: functions initialize and delete simulation fields
  * **init_del_ghost_container.cpp**: functions initialize and delete interface ghost data.
  * **init_del_fd_coeff.cpp**: functions initialize and delete finite difference coefficients.
  * **init_src.cpp**: functions initialize a Ricker source.
  * **init_del_rec_trace**: functions initialize and delete a set of receiver traces.
  * **init_del_movie.cpp**: functions initlalize and delete rsf movies.
  * **init_del_res_trace.cpp**: functions initialize and delete residue traces
  * **update_*.cpp**: functions update velocity/stresses/ghost data at grid interface.
  * **compute_energy.cpp/compute_res.cpp**: functions compute energy and residue at a composite grid.
  * **print.cpp**: functions print receiver data/residues/movies to rsf files.
* **test_utils.hpp**: list of test function signatures with functions defined in:
  * **test_sbp.cpp**: functions test SBP property.
* **tests**: test scripts

===================================================================

## Future Work
- [ ] Extend this composite grid scheme to deal with fluid-solid interface.
- [ ] Extend this composite grid scheme to high order FDTD.
- [ ] MPI parallelization.

