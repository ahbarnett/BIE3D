# BIE3D: MATLAB tools for boundary integral equations on surfaces in 3D

This is a preliminary set of high-order accurate quad-panel based surface quadratures for kernels that have on-surface weak singularities no more singular than 1/r. It currently contains only a torus and its modulation via a general smooth height function, with uniform arbitrary-order quad patches, with the Laplace (an elliptic BVP) and wave-equation (hyperbolic BVP) kernels.

Main author:  Alex Barnett

Contributions: Tom Hagstrom - f90 modules for interpolation from time grid.

Version: 20181110

### Dependencies

MATLAB. Codes have not been tested on MATLAB versions prior to R2012a.

For `timedomainwaveeqn`:

* Fortran compiler to build Hagstrom time interpolation and MEX interface.  
* Optionally: `fsparse` from [stenglib](https://github.com/stefanengblom/stenglib), compiled with `make('openmp',true)`, for fast multithreaded sparse matrix assembly.  

### Installation

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE3D`) directory, and run `bie3dsetup` to add all needed directories to your path. 

Test by running `testall` which tests the Laplace case on a torus, and should produce lots of error outputs close to machine precision, convergent sequences of numbers, and some plots, and yet not crash.


### Directories

`kernels`  : Laplace evaluation, including on-surface (self-eval)  
`surfaces` : smooth surface generators  
`singquad` : special surface quadratures for weakly singular kernels  
`utils`    : general numerical and plot utilities  
`test`     : test codes (other than built-in self-tests)  
`timedomainwaveeqn` : TDBIE for wave equation, in progress  
`solvers` : 2D BVP solver example codes, also serve to test kernels  
`doublyperiodic` : codes for doubly-periodic geometries  
