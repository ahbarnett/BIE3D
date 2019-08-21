function [h h2] = showsurf(s,c,o)
% SHOWSURF.  Plot nodes in a 3D surface struct, and possibly normals.
%
% h = showsurf(s) plots in 3D, on current figure, the points s.x,
%  with normals given by s.nx, returning a graphics handle h to the points.
%
% h = showsurf(s,c) uses the color character in c.
%
% h = showsurf(s,c,o) controls various options:
%  o.normals = 0 (no normals), 1 (show normals, default).
%  o.alpha in [0,1], adds alpha-transparent surface (alpha=1 opaque), in the
%   case of torus-like topology. [h h2] = showsurf(...) returns h2 surf handle.
%
% See for test: self-test for setup_torus_doubleptr

% Barnett 8/21/19
if nargin<2, c = 'k'; end
if nargin<3, o=[]; end
if ~isfield(o,'normals'),o.normals=1; end  

h = plot3(s.x(1,:),s.x(2,:),s.x(3,:),[c '.'],'markersize',2);       % pts
if o.normals
  hold on;
  y = s.x + 0.05*s.nx;
  plot3([s.x(1,:);y(1,:)],[s.x(2,:);y(2,:)],[s.x(3,:);y(3,:)],[c '-']);
end  
if isfield(o,'alpha')    % add surface
  [X,Y,Z] = deal(nan(s.Na+1, s.Nb+1));
  X(1:s.Na,1:s.Nb) = reshape(s.x(1,:),[s.Na s.Nb]);
  Y(1:s.Na,1:s.Nb) = reshape(s.x(2,:),[s.Na s.Nb]);
  Z(1:s.Na,1:s.Nb) = reshape(s.x(3,:),[s.Na s.Nb]);  
  X(:,end) = X(:,1); X(end,:) = X(1,:); X(end,end) = X(1,1);
  Y(:,end) = Y(:,1); Y(end,:) = Y(1,:); Y(end,end) = Y(1,1);
  Z(:,end) = Z(:,1); Z(end,:) = Z(1,:); Z(end,end) = Z(1,1);
  linecol = get(h,'Color'); linecol = linecol(:)';  % row vec
  C = kron(linecol,ones(s.Na+1,s.Nb+1));   % same color all vertices
  C = reshape(C,[s.Na+1,s.Nb+1,3]);
  h2 = surf(X,Y,Z,C,'FaceAlpha',o.alpha,'edgecolor','none');
  set(h2,'ambientstrength',0.7,'facelighting','gouraud');
end  
axis equal vis3d
xlabel('x'); ylabel('y'); zlabel('z');
