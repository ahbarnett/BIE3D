# BIE3D: MATLAB tools for boundary integral equations on surfaces in 3D

This contains some research codes for high-order accurate 3D Nystrom BIE for scalar elliptic PDEs with constant coefficients inside or outside of a surface.
In more detail, it has: global double periodic trapezoid rule and quad-panel based surface quadratures for kernels that have on-surface weak singularities no more singular than 1/r. For the torus and its modulation via a general smooth radius function it has on-surface quadratures only for uniform arbitrary-order quad patches, for the Laplace (an elliptic BVP) and wave-equation (hyperbolic BVP) kernels. For smooth deformations of the sphere, it has global QFS quadratures for the Laplace kernel, built on global spectral interpolations. It also hosts an old self-contained doubly-periodic Laplace dipole summation code.

Author:  Alex Barnett

Contributions: Tom Hagstrom - f90 modules for interpolation from time grid.

Version: 20230317 (tested on MATLAB R2022a)

![cruller with panel quadrature and on-surface singular scheme](pics/cruller_panelquad.png)

### Dependencies

MATLAB. Codes have not been tested on MATLAB versions prior to R2012a.

For `timedomainwaveeqn`:

* Fortran compiler to build Hagstrom time interpolation and MEX interface.
* Some driver scripts need you to have the MATLAB/Octave tool [memorygraph](https://github.com/ahbarnett/memorygraph)  
* Optionally: `fsparse` from [stenglib](https://github.com/stefanengblom/stenglib), compiled with `make('openmp',true)`, for fast multithreaded sparse matrix assembly.  

### Installation

Download using `git` or as a zip (see green button above).

Open MATLAB in the top level (`BIE3D`) directory, and run `bie3dsetup` to add all needed directories to your path. 

Test by running `testall` which currently tests Laplace quadratures on a torus, and should produce lots of error outputs close to machine precision, convergent sequences of numbers, and some plots, taking around 30-60 secs total.

### Directories

`kernels`  : Laplace evaluation, including on-surface (self-eval)  
`surfaces` : smooth surface generators  
`singquad` : special surface quadratures for weakly singular kernels  
`utils`    : general numerical and plot utilities  
`test`     : test codes (other than built-in self-tests)  
`timedomainwaveeqn` : time-domain integral-equations for acoustics codes (see [movie](http://users.flatironinstitute.org/~ahb/images/cruller_scatt_plane_pulse_m4_p6_np24_hi.mp4))  
`doublyperiodic` : an old self-contained code for doubly-periodic Laplace dipoles in 3D  
