function [S D Dp] = tdSDmats(x,y,ny,w)
% TDSDMATS   fill rows applying time-domain S,D kernels (w/ wei) to ret dens
%
% [S D Dp] = tdSDmats(x,y,ny,w)
%
% Inputs:
%  x - targets (3-by-M)
%  y,ny - sources and source normals (both 3-by-N)
%  w - source quadrature weights (1-by-N)
% Outputs:
%  S,D,Dp - MxN matrices which, when each row is dotted into a vector of
%           corresponding retarded densities,
%           computes the quadrature approx to the potential at each target x.
%           ie, SLP = s*sigma(y), and DLP = d*tau(y) + dp*tau_t(y).
%
% Notes:
% 1) this duplicates but combines the matrix-filling codes from
%    {Lap,Ret}{S,D}eval_panels.
% 2) not tied to any particular quadr scheme; can also be used for aux nodes.
%    TODO: make accept panel arrays t, s ?
% 3) dense, so don't call with too large a product NM.
% 4) meaningless to use for plain matvecs since the retardation in the dens vec
%    would depend on the target point x.
%
% See also: TEST_TDGRF.m which tests it (choose side=+-1)

% Barnett 12/29/16
M = size(x,2);
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
rr = d1.^2+d2.^2+d3.^2;              % M*N
ir = 1./sqrt(rr);
S = ir .* repmat((1/4/pi)*w(:)', [M 1]);    % right-mult w/ wei & 1/4pi fac
ddotn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
Dp = ddotn.*ir .* S;       % already have correct wei and prefacs
D = ir.*Dp;
