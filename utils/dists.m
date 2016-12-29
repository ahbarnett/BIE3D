function r = dists(x,y)
% DISTS  Euclidean distances from one to many points in R^d.
%
% r = dists(x,y)
%  x : dx1 target
%  y : dxN sources
%  r : 1xN row vec of distances ||x-y_j||,  j=1,..,N
% Note: alternatively, x and y can be swapped. But at least one of x,y must be
% size dx1 (ie, a single point).

% Barnett 12/29/16
r = sqrt(sum(bsxfun(@minus,x,y).^2,1));
