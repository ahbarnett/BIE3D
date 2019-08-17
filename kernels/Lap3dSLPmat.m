function [A An] = Lap3dSLPmat(t,s)
% LAP3DSLPMAT.  dense Laplace SLP Nystrom eval matrix from sources to targets
%
% [A An] = Lap3dSLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% NB: fills M*N matrices, so assumes small cases. Faster if An not needed.
%
% Without args, does partial math test.
%
% For full test, see: ../test/test_LapGRF_torus_global.m
%
% Todo: * replace with C+omp MEX, combine with DLP case if needed together.

% Barnett 8/16/19

if nargin==0, test_Lap3dSLPmat; return; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
A = bsxfun(@times, 1./sqrt(rr), s.w*(1/4/pi));  % including src quadr wei
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  An = bsxfun(@times, ddottn./(sqrt(rr).*rr), s.w*(1/4/pi));  % monopole deriv, incl src quad wei
end

function test_Lap3dSLPmat
% math test for grad only: Gauss's Law on torus surface...
t = setup_torus_doubleptr(1,0.5);
xin = [0.9; -0.2; 0.1]; xout = [1.9; 0.7; 1.0];    % both "far" from surf
s.x = [xin,xout]; s.w = 1;                         % source pointset
[~,An] = Lap3dSLPmat(t,s);                         % this test ignores val
fluxes = t.w * An;                                 % w is row, leaves row 2-vec
fprintf('torus SLP An errs: int %.3g, ext %.3g\n',fluxes(1)-1,fluxes(2))
