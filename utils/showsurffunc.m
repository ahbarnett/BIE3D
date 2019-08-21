function [h h2] = showsurffunc(s, f, o)
% SHOWSURFFUNC   show real-valued scalar func on surface (torus-like for now)
%
% [h h2] = showsurffunc(s, f, o) plots a 3D surface with colorscale indicating
%  value of the function f defined on the nodes of the quadrature in s.
%
% Inputs:
%  s = either: cell array of quad panels (panel case; only torus-like surface
%              for now)
%      or: global surface struct, with fields x, etc.
%  f = list (any shape) of real function values on the nodes s.x.
%  Automatically decides if 1-sided or symmetrized colorscale is best.
% Outputs: h - handle to surface, h2 - handle to colorbar
%
% o controls opts:
%  o.sc = sets color scale
%  o.nofig = if is a field, write on current fig

% Barnett 7/19/16, based on 2/21/13 qp3d (swapped np,mp, killed 2d plot)
% added global case 8/21/19
if nargin==0, test_showsurffunc; return; end
if nargin<3, o = []; end
if ~isfield(o,'nofig'), newfig = 1; else newfig = 0; end

f = real(f);
if ~isfield(o,'sc'), o.sc = max(abs(f)); end % colorscale, then symm or not...
if sum(f<0)==0, cs = [0 o.sc]; else cs = o.sc*[-1 1]; end

% set up arrays, including wrapping to close the surface...
e = 1;   % 1 if extra row & co needed for sealing up (closed surf)
if isstruct(s)                           % global, torus
  fr=reshape(f,[s.Na s.Nb]);
  [X,Y,Z,f] = deal(nan(s.Na+e, s.Nb+e));  % NB f gets overwritten
  size(s.x(1,:))
  size(X(1:s.Na,1:s.Nb))
  X(1:s.Na,1:s.Nb) = reshape(s.x(1,:),[s.Na s.Nb]);
  Y(1:s.Na,1:s.Nb) = reshape(s.x(2,:),[s.Na s.Nb]);
  Z(1:s.Na,1:s.Nb) = reshape(s.x(3,:),[s.Na s.Nb]);
elseif strcmp(s{1}.topo,'torus')         % panels, the case we know: torus
  np = s{1}.np; mp = s{1}.mp; q=s{1}.p; x = s{1}.t(2,1:q); n=size(s{1}.t,2);
  npan = mp*np;
  panind = reshape(1:npan,[mp np]); % 2d grid of panel numbers (sets dof order)
  fr = nan(q*mp, q*np); % reorder f vector to a 2d global grid fr...
  for j=1:np, for i=1:mp, k=panind(i,j);   % panel number
      fr((i-1)*q+(1:q),(j-1)*q+(1:q)) = reshape(f(s{k}.inds),q,q);
  end, end
  % seal it up so looks periodic
  [X,Y,Z,f] = deal(nan(q*mp+e, q*np+e));  % NB f gets overwritten
  for j=1:np, for i=1:mp, k=panind(i,j);   % panel number
      X(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(1,:),q,q);
      Y(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(2,:),q,q);
      Z(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(3,:),q,q);
    end, end
end
if e  % add in repeated edges to seal up
  f(1:end-1,1:end-1) = fr;
  X(:,end) = X(:,1); X(end,:) = X(1,:); X(end,end) = X(1,1);
  Y(:,end) = Y(:,1); Y(end,:) = Y(1,:); Y(end,end) = Y(1,1);
  Z(:,end) = Z(:,1); Z(end,:) = Z(1,:); Z(end,end) = Z(1,1);
  f(:,end) = f(:,1); f(end,:) = f(1,:); f(end,end) = f(1,1); 
else  % interface ? - TODO make seal up poss w/ Bloch phases
  f = fr;
end
% plot...
if newfig, figure; end   %else, hold on; end
h = surf(X,Y,Z, f); set(h,'ambientstrength',0.7,'facelighting','gouraud')
shading interp; xlabel('x'); ylabel('y'); zlabel('z');
if newfig, lightangle(45,0); end
axis equal vis3d; 
h2 = colorbar; colormap(jet(256)); caxis(cs);

%%%%%%%%%%%%
function test_showsurffunc
% panel case...
type = 'torus'; so.a = 1; so.b = 0.5; [s N] = create_panels(type,so);
x = getallnodes(s);
f = sin(1.7*x(1,:)-0.5*x(2,:)+0.9*x(3,:))';    % col vec of func vals
showsurffunc(s, f); title('panel torus');
% global case...
s = setup_torus_doubleptr(so.a,so.b);
f = sin(1.7*s.x(1,:)-0.5*s.x(2,:)+0.9*s.x(3,:))';
showsurffunc(s, f); title('double-PTR torus');
