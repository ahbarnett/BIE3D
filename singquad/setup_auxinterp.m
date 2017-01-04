function I = setup_auxinterp(tt,o)
% SETUP_AUXINTERP   build interpolation matrices from near panels to aux nodes
%
% I = setup_auxinterp(tt,o) returns a struct giving interp matrices from smooth
%     quadrature to auxiliary nodes, all in standard panels [-1,1] in param
%     space.
% Assumes tt is tensor product nodes in param space.
%
% Inputs:
%   tt - (2-by-p^2) list of param coords for the target pts. The targets are
%        assumed same for each of the panels in the 3x3 block.
%   o - opts to pass to panel_sing_auxquad
%
% Source panel ordering  2 5 7      since the self comes first, then neighbors
%                        3 1 8      (note this shows the parameter plane)
%                        4 6 9
%
% Outputs: I.L - 9-by-1 cell array of interp matrices from n=p^2 nodes in
%                  each near src panel to all the aux nodes in that panel
%          I.auxinds - 9-by-1 cell array of index lists, each giving which aux
%                  nodes (indexed from the full list, which is in naux-fast,
%                  targ-node-slow ordering) are in this source panel.
%          I.auxindsbytarg - 9-by-p^2 cell array of index lists for which
%                  aux nodes indices are in each src panel, for each targ node.
%          I.Lbytarg - 9-by-p^2 cell array of interp matrices from p^2 nodes in
%                  each src panel to aux nodes in that panel, for each of p^2
%                  targets (the 2nd cell array dim indexes target).
%          I.z  - (2-by-(naux*p^2)) the aux nodes used, in param space,
%                 in the ordering used 
%          I.cen - 9-by-1 cell array of 2-by-1 centers of the source boxes
%                 in param space.
%
% See also: CREATE_PANELS, where the std nei ordering must match.

% Barnett 7/17/16, bytarg needed by tdSDinterp... 1/3/17
if nargin==0, test_setup_auxinterp; return; end

z = panel_sing_auxquad(tt,o);   % slow but ok to do it again
n = size(tt,2);            % smooth nodes per panel, = p^2
p = sqrt(n); if p~=round(p), error('n is not p^2 for some integer p!'); end
t1 = tt(1,1:p:1+(p-1)*p); t2 = tt(2,1:p);    % 1d nodes (should be gauss!)
naux = size(z,2);           % aux nodes per target pt
I.z = reshape(z, [2 naux*n]);   % horiz stack all aux node locs in param space
clear z
I.L = cell(9,1); I.auxinds = cell(9,1); I.cen = cell(9,1);

% x and y centers of 9 cells in std ordering...
cen = 2*[0 -1 -1 -1 0 0 1 1 1; 0 -1 0 1 -1 1 -1 0 1];  % incr in y fast, then x
for i=1:9          % loop over self and neighbors
  xc = cen(1,i); yc = cen(2,i);  % this box center x y
  x = bsxfun(@minus, I.z, [xc;yc]);  % all aux nodes rel to box ctr
  I.auxinds{i} = find(abs(x(1,:))<1 & abs(x(2,:))<1);  % keep those in [-1,1]^2
  x = x(:,I.auxinds{i});
  tind = kron(1:n,ones(1,naux)); tind = tind(I.auxinds{i}); % target indices
  I.L{i} = tensorprod_interp(x,t1,t2);
  I.cen{i} = cen(:,i);
  for j=1:n                     % loop over targ nodes
    I.auxindsbytarg{i}{j} = I.auxinds{i}(tind==j);
    I.Lbytarg{i}{j} = I.L{i}(tind==j,:);  % maps p^2 src to aux_ij
  end
end

%%%%%%
function test_setup_auxinterp
p = 8;   % build a single panel
x = gauss(p);   % Gauss-Legendre nodes & weights on [-1,1]
[x1 x2] = meshgrid(x); tt = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
o.nt = 16; o.nr = 8;   % # aux nodes in theta,rho
tic
I = setup_auxinterp(tt,o);  % slow!
toc
f = @(x) sin(x(1,:) + x(2,:)/2 + 0.7).*exp(-0.2*x(1,:) + 0.6*x(2,:));  % smooth in [-1,1]^2, maps cols in R2 to vals (row output)
data = f(tt)';    % col vec of p^2 samples (smooth quadr). Same f for each i:
for i=1:9, L = I.L{i}; ii = I.auxinds{i};
  z = bsxfun(@minus, I.z(:,ii), I.cen{i});
  err = max(abs(L*data - f(z)'));
  fprintf('near cell %d: max interp err to aux nodes in that cell = %.3g\n',i,err)
end
