function [d x y i j] = dist_pointclouds(b1,b2,meth)
% DIST_POINTCLOUDS find min dist between two sets of points
%
% [d x y i j] = dist_pointclouds(b1,b2) where b1, b2 are objects with
%  fields x giving 3xN array of points in R3, returns
%  d = dist(b1,b2) an estimate of the closest distance, and x,y the pair
%  achieving this distance, and i, j their respective column
%  indices in b1.x, b2.x.
%
% [d x y i j] = dist_pointclouds(b1,b2,meth) controls method:
%   meth='b' brute-force O(N^2), guaranteed
%   meth='g' greedy O(N), often gets stuck in poor minima (default)
%
%  Note: for surfaces it is only as accurate as the discretization!
%
% Called without arguments it does a self-test (see code for example usage).
%
% See also: dists

% Barnett 8/15/21
if nargin==0, test_dist_pointclouds; return, end
if nargin<3, meth='b'; end

N1 = size(b1.x,2);

if meth=='b'
  if N1*size(b2.x,2)>1e9, warning('brute force >1e9 elements!'); end
  dd = (b1.x(1,:) - b2.x(1,:)').^2 + (b1.x(2,:) - b2.x(2,:)').^2 + (b1.x(3,:) - b2.x(3,:)').^2;
  [d,ii] = min(dd(:));                     % note d = dist^2
  j = mod(ii,N1); i = floor(ii/N1);
  x = b1.x(:,i); y = b2.x(:,j);
  
elseif meth=='i'
  i=randi(N1);
  x = mean(b1.x,2);    % center of one
  maxits = 20;
  for k=1:maxits
    xold = x;
    [d j] = min(sum((b2.x - x).^2,1));         % note d = dist^2
    y = b2.x(:,j);
    [d i] = min(sum((b1.x - y).^2,1));
    x = b1.x(:,i);
    if norm(x-xold)==0, break, end
  end
  
else
  error('unknown meth!');
end
d = sqrt(d);


%%%%%%%%%%%
function test_dist_pointclouds
res = [100 50];
%res = 2.5*[100 50];

b = ellipsoid(1,2,3); b = setupsurfquad(b,res);
b1.x = b.x;   % make sure just this field used

Rz = @(t) [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];  % z-ax rot mat
Ry = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];  % y-ax rot mat
rng(0);
R = Rz(2*pi*rand) * Ry(pi*rand) * Rz(2*pi*rand);    % Euler angles
b2.x = R*b.x + [2;2;-1];              % rotation then translate

tic; [d x y i j] = dist_pointclouds(b1,b2,'b'); toc  % brute force
d

z=[b1.x,b2.x]; figure; plot3(z(1,:),z(2,:),z(3,:),'k.'); axis vis3d equal
% line connecting nearest pts...
hold on; plot3([x(1) y(1)],[x(2) y(2)],[x(3) y(3)], '.-', 'markersize',20);
view(40,40);
