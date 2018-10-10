function [jmax jmin umat upmat] = interpmat(tinterp,dt,m)
% INTERPMAT  return vectors of coeffs which interp from reg t grid to target t's
%
% [jmax jmin umat upmat] = interpmat(tinterp,dt,m)
%
% Inputs:
% tinterp - list of desired times, assumed negative (the current time = 0)
% dt      - time step
% m       - order of accuracy: degree 2m+1 spline
%
% Outputs:
% jmax    - maximum time index for interpolation data - <= 0 with current time 0
% jmin    - minimum time index for interpolation data - <= 0 with current time 0
% umat    - matrix of dimension r * (jmax-jmin+1) for u interpolation, r is number of
%           times in the tinterp input list.
% upmat   - matrix of dimension r * (jmax-jmin+1) for du/dt interpolation
%
% Notes: MEX interface to code by Tom Hagstrom, July 2016.
% See TestInterpMat executable, made by make exec.
% Interface: Alex Barnett 12/19/16

if m==0, error('m must be positive!'); end
tinterp = tinterp(:); r = numel(tinterp);
if sum(tinterp>0), warning('some tinterps are positive!'); end

% preallocate outputs via upper bnd on nc, number of output indices...
maxinds = 5*m + ceil((max(tinterp)-min(tinterp))/dt);
umat = zeros(r,maxinds); upmat = umat;     % allocate outputs for fortran

% note this is to a simpler non-dynamically allocated umat,upmat version...
mex_id_ = 'interpmatnoalloc(i int, i double[], i double, i int, o int[x], o int[x], io double[], io double[], i int)';
[jmax, jmin, umat, upmat] = gateway(mex_id_, r, tinterp, dt, m, umat, upmat, maxinds, 1, 1);

nc = jmax-jmin+1;     % desired number of output grid pts
umat = umat(:,1:nc);  % resize to correct size...
upmat = upmat(:,1:nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
