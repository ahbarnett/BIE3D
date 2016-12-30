function r = dists(x,y)
% DISTS  Euclidean distance matrix from many to many points in R^3.
%
% r = dists(x,y)
%  x : 3xM target
%  y : 3xN sources
%  r : MxN matrix of distances r_{ij} = ||x_i-y_j||

% Barnett 12/29/16
%r = sqrt(sum(bsxfun(@minus,x,y).^2,1)); % only 1 src or 1 targ, but R^d.
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
r = sqrt(d1.^2+d2.^2+d3.^2);          % matrix, size M*N
