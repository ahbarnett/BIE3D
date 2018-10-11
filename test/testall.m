% complete set of tests for BIE3D. All should produce small numbers &
% occasional figures and tables.
% Barnett 7/19/16

% utils
checkgrad
interpmat_1d
tensorprod_interp

% surfaces
create_panels
test_wobblytorus

% singquad
panel_sing_auxquad
add_panels_auxquad
setup_auxinterp

% kernels
LapDeval_panels
LapSeval_panels

% test
test_LapGRF

% solvers

% doublyperiodic

% timedomainwaveeqn
