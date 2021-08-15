function [d x y its] = dist_ellipsoids(E1,R1,t1,E2,R2,t2, abstol)
% DIST_ELLIPSOIDS  find min dist between two ellipsoids
%
% [d x y] = dist_ellipsoids(E1,R1,t1,E2,R2,t2)
%  finds dist between two ellipsoids, each defined by semiaxis list E,
%  3x3 rotation matrix R, and center vector t.
%  d an estimate of the closest distance, and x,y the pair achieving it
%  (x in 1st ellipsoid, y in 2nd).
%
% [d x y its] = dist_ellipsoids(E1,R1,t1,E2,R2,t2, abstol)
% also controls estim absolute tolerance for dist
% (crappy, need to set to tol^2 !)
% and reports # iters.

% Called without arguments it does a self-test (see code for example usage).

% Barnett 8/15/21
if nargin==0, test_dist_ellipsoids; return, end
if nargin<7, abstol = 1e-6; end

t = t2(:)-t1(:);   % col vec
R = inv(R1)*R2;
iR = inv(R);
% now E1 is oriented vertically at origin, E2 is at t with rot mat R.

x = t;    % center of E2
dold = inf;
maxit = 1e4;          % takes ages since slow decay when close!
for k=1:maxit
  [d y] = dist_ellipsoid(x,E1);
  [d x] = dist_ellipsoid(iR*(y-t),E2);   % solve rel to coords at E2
  x = t+R*x;            % convert back to E1 coord sys
  if abs(d-dold)<abstol, break, end
  dold = d;
%  k,x,d
end
its = k;

%%%%%%%%%%%
function test_dist_ellipsoids

E1 = [1 2 3];   % ellipsoid 1, semiaxes
R1 = eye(3);
t1 = zeros(3,1);

E2 = [1 2 3];       % ellipsoid 2
Rz = @(t) [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];  % z-ax rot mat
Ry = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];  % y-ax rot mat
rng(0);
R2 = Rz(2*pi*rand) * Ry(pi*rand) * Rz(2*pi*rand);    % Euler angles
t2 = [2;2;-1];              % translate

tic; [d0 x0 y0 its] = dist_ellipsoids(E1,R1,t1,E2,R2,t2, 1e-10);  % "exact"
d0, its, toc
%norm(x0-y0)

res = [100 50];     % compare against crappier discr method...
b = ellipsoid(E1(1),E1(2),E1(3)); b=setupsurfquad(b,res);
b1.x = R1*b.x + t1;
b = ellipsoid(E2(1),E2(2),E2(3)); b=setupsurfquad(b,res);
b2.x = R2*b.x + t2;

tic; [d x y i j] = dist_pointclouds(b1,b2,'b');  % brute force
d, toc

z=[b1.x,b2.x]; figure; plot3(z(1,:),z(2,:),z(3,:),'k.'); axis vis3d equal
% line connecting nearest pts... red is exact, blue the crappy discr ans:
hold on; plot3([x(1) y(1)],[x(2) y(2)],[x(3) y(3)], '.-', 'markersize',20);
plot3([x0(1) y0(1)],[x0(2) y0(2)],[x0(3) y0(3)], 'r.-', 'markersize',20);

view(40,0);  % zoom in on it
set(gca,'CameraTarget',x0);
set(gca,'Position',[0.130000000000000   0.110000000000000   0.775000000000000   0.815000000000000]);
