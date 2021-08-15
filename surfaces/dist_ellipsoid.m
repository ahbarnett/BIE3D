function [dist,z,lambda] = dist_ellipsoid(y,semiaxes)
% DIST_ELLIPSOID  distance from point in 3D to standard ellipsoid
%
% [dist,z,lambda] = dist_ellipsoid(x,semiaxes) where semiaxes=[a b c], and
%  x is point in R3,
% returns dist, point achieving it z, and lagrange multiplier lambda
% (lambda<0 if was inside)
%
% Without arguments does self-test.

% Code tweaked by Barnett 8/15/21, from Beni Bogosel's of 2013.
% https://mathproblems123.wordpress.com/2013/10/17/distance-from-a-point-to-an-ellipsoid/
if nargin==0, test_dist_ellipsoid; return; end

y = y(:);
vec = 1./semiaxes.^2;
vec = vec(:);

epst = 1e-12;
x0=0;
x1=0.1;  % why?
while abs(x1-x0)>epst
x0 = x1;
x1 = x0-fun(x0,y,vec)/der(x0,y,vec);    % Newton
end
lambda = x1;
%z = (eye(length(vec))+lambda*diag(vec))\y;   % lin solve is diagonal
z = y ./ (1+lambda*vec);
dist = norm(y-z);
if lambda<0, dist=0; end         % handle inside case as zero-dist

function res = fun(lam,y,vec);
res = sum(vec.*(y./(lam*vec+1)).^2)-1;
 
function res2 = der(lam,y,vec);
res2 = -2*sum((vec.^2).*(y.^2)./((lam*vec+1).^3));

%%%%%%%%%%
function test_dist_ellipsoid
x = [1 2 3];
d=dist_ellipsoid(x,[1 1 1]);  % sphere case
d - (norm(x)-1)

[d,z,lam] =dist_ellipsoid([0 0 5],[1 2 3]); d-2   % axes cases
d=dist_ellipsoid([0 4 0],[1 2 3]); d-2
d=dist_ellipsoid([3 0 0],[1 2 3]); d-2

[d,z,lam] =dist_ellipsoid([0 0 2],[1 2 3]); d  % inside, should give 0
[d,z,lam] =dist_ellipsoid([0 0 -2],[1 2 3]); d  % "
