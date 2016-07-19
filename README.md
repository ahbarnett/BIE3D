# BIE3D: MATLAB tools for boundary integral equations on surfaces in 3D

This is a preliminary set of high-order accurate quad-panel based surface quadratures for kernels that have on-surface weak singularities no more singular than 1/r. It currently contains only a torus, and the Laplace kernels.

It is intended to include a sandbox including the time-domain BIE for the wave equation.

Main author:  Alex Barnett

Contributions: Tom Hagstrom - interpolation from time grid.

Version: 20160719

### Installation

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE3D`) directory, and run `bie3dsetup` to add al
l needed directories to your path. 

Test by running `testall` which should produce lots of error outputs close to ma
chine precision, convergent sequences of numbers, and some figures, and yet not crash.

Codes have not been tested on MATLAB versions prior to R2012a.

### Directories

`kernels`  : Laplace evaluation, including on-surface (self-eval)  
`surfaces` : smooth surface generators  
`singquad` : special surface quadratures for weakly singular kernels  
`utils`    : general numerical and plot utilities  
`test`     : test codes (other than built-in self-tests)  
`timedomainwaveeqn` : TDBIE for wave equation (not yet)  
`solvers` : 2D BVP solver example codes, also serve to test kernels  
`doublyperiodic` : codes for doubly-periodic geometries  
