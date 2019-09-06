% complete set of tests for BIE3D. All should produce small numbers &
% occasional figures and tables.
% Barnett 7/19/16 for panel stuff; updated 8/22/19 for global quad stuff.

% utils (directory)
checkgrad
interpmat_1d
tensorprod_interp
peri2dspecinterp
intxperiinterp
intxperiinterpmat
showsurf
showsurffunc
setupdoubleptr
setupspherequad

% surfaces
setup_torus_doubleptr
create_panels
test_wobblytorus

% singquad
panel_sing_auxquad
add_panels_auxquad
setup_auxinterp

% kernels
Lap3dSLPmat
Lap3dDLPmat
LapDeval_panels
LapSeval_panels

% test
test_LapGRF_torus_global
test_LapGRF_torus_panels

% solvers

% doublyperiodic

% timedomainwaveeqn
