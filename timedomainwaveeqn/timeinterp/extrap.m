function xcof = extrap(m)
% EXTRAP  row vector of coeffs to extrapolate one pt beyond regular grid
%
% xcof = extrap(m)
%
% Inputs:
%  m - order
% Outputs:
%  xcof - row vector of m+1 real coefficients to be applied to regular grid
%         to extrapolate one grid point beyond the end (ie to m+1, if the
%         existing grid is 0,1,..,m)
%
% Notes: MEX interface to code by Tom Hagstrom, July 2016.
% Interface: Alex Barnett 12/23/16

xcof = zeros(1,m+1);
mex_id_ = 'extrapnoalloc(io double[], i int)';
[xcof] = gateway(mex_id_, xcof, m);
