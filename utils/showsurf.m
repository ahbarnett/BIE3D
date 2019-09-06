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
%  o.alpha in [0,1], adds alpha-transparent surface (alpha=1 opaque), if it
%  can read s.topo giving the topology (torus, sphere,...)
% [h h2] = showsurf(...) returns h2 surf handle.
%
% For test: run without arguments.

% Barnett 8/21/19, 9/3/19 for sphere topo.
if nargin<1, test_showsurf; return; end
if nargin<2, c = 'k'; end
if nargin<3, o=[]; end
if ~isfield(o,'normals'),o.normals=1; end  

h = plot3(s.x(1,:),s.x(2,:),s.x(3,:),[c '.'],'markersize',2);       % pts
if o.normals & isfield(s,'nx')
  hold on;
  y = s.x + 0.1*s.nx;
  plot3([s.x(1,:);y(1,:)],[s.x(2,:);y(2,:)],[s.x(3,:);y(3,:)],[c '-']);
end  
axis equal vis3d
xlabel('x'); ylabel('y'); zlabel('z');
set(gca,'clipping','off');
if isfield(o,'alpha')    % add surface
  if s.topo=='t'         % torus-like
    nu = s.Nu+1; nv = s.Nv+1;
    [X,Y,Z] = deal(nan(nu,nv));
    X(1:s.Nu,1:s.Nv) = reshape(s.x(1,:),[s.Nu s.Nv]);
    Y(1:s.Nu,1:s.Nv) = reshape(s.x(2,:),[s.Nu s.Nv]);
    Z(1:s.Nu,1:s.Nv) = reshape(s.x(3,:),[s.Nu s.Nv]);  
    X(:,end) = X(:,1); X(end,:) = X(1,:); X(end,end) = X(1,1);
    Y(:,end) = Y(:,1); Y(end,:) = Y(1,:); Y(end,end) = Y(1,1);
    Z(:,end) = Z(:,1); Z(end,:) = Z(1,:); Z(end,end) = Z(1,1);
  elseif s.topo=='s'     % sphere-like
    if numel(s.Nu)==1    % tensor-prod case
      nu = s.Nu+1; nv = s.Nv+2;
      [X,Y,Z] = deal(nan(nu,nv));
      X(1:nu-1,2:nv-1) = reshape(s.x(1,:),[s.Nu s.Nv]);
      Y(1:nu-1,2:nv-1) = reshape(s.x(2,:),[s.Nu s.Nv]);
      Z(1:nu-1,2:nv-1) = reshape(s.x(3,:),[s.Nu s.Nv]);
      npole = s.Z(0,1); spole = s.Z(0,-1);     % assumes analytic avail
      X(:,1) = spole(1); Y(:,1) = spole(2); Z(:,1) = spole(3);   % patch the poles
      X(:,end) = npole(1); Y(:,end) = npole(2); Z(:,end) = npole(3);
      X(end,:) = X(1,:);   % wrap in u
      Y(end,:) = Y(1,:);
      Z(end,:) = Z(1,:);
    else
      %disp('(dont know how to show non-tensor-prod sphere-like as surf)')
      return
    end
  end
  linecol = get(h,'Color'); linecol = linecol(:)';  % row vec
  C = kron(linecol,ones(nu,nv));   % same color all vertices
  C = reshape(C,[nu,nv,3]);
  h2 = surf(X,Y,Z,C,'FaceAlpha',o.alpha,'edgecolor','none');
  set(h2,'ambientstrength',0.7,'facelighting','gouraud');
end  

%%%%%%
function test_showsurf   % exploit other's self-tests
setupdoubleptr
setupspherequad
