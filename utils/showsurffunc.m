function h = showsurffunc(s, f, o)
% SHOWSURFFUNC   show real-valued vector on surface (torus-like for now)
%
% h = showsurffunc(s, f, sc, o)
% Inputs:
%  s = cell array of quad panels (only torus-like surface for now)
%  f = list of real function values on the nodes s.x.
%  Automatically decides if 1-sided or symmetrized colorscale is best.
%
% o controls opts:
%  o.sc = sets color scale
%  o.nofig = if a field, write on current fig

% Barnett 7/19/16, based on 2/21/13 qp3d (swapped np,mp, killed 2d plot)
if nargin==0, test_showsurffunc; return; end
if nargin<3, o = []; end
if ~isfield(o,'nofig'), fig = 1; else fig = 0; end

f = real(f);
if ~isfield(o,'sc'), o.sc = max(abs(f)); end % colorscale, then symm or not...
if sum(f<0)==0, cs = [0 o.sc]; else cs = o.sc*[-1 1]; end

if strcmp(s{1}.topo,'torus')             % the only case we know
  np = s{1}.np; mp = s{1}.mp; q=s{1}.p; x = s{1}.t(2,1:q); n=size(s{1}.t,2);
  npan = mp*np;
  panind = reshape(1:npan,[mp np]); % 2d grid of panel numbers (sets dof order)
  fr = nan(q*mp, q*np); % reorder f vector to a 2d global grid fr...
  for j=1:np, for i=1:mp, k=panind(i,j);   % panel number
      fr((i-1)*q+(1:q),(j-1)*q+(1:q)) = reshape(f(s{k}.inds),q,q);
  end, end
  % seal it up so looks periodic
  e = 1;   % 1 if extra row & co needed for sealing up.
  X = nan(q*mp+e, q*np+e); Y = X; Z = X; f = X;
  for j=1:np, for i=1:mp, k=panind(i,j);   % panel number
      X(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(1,:),q,q);
      Y(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(2,:),q,q);
      Z(q*(i-1)+(1:q),q*(j-1)+(1:q))=reshape(s{k}.x(3,:),q,q);
  end, end
  if e % add in repeated edges to seal up
    f(1:end-1,1:end-1) = fr;
    X(:,end) = X(:,1); X(end,:) = X(1,:); X(end,end) = X(1,1);
    Y(:,end) = Y(:,1); Y(end,:) = Y(1,:); Y(end,end) = Y(1,1);
    Z(:,end) = Z(:,1); Z(end,:) = Z(1,:); Z(end,end) = Z(1,1);
    f(:,end) = f(:,1); f(end,:) = f(1,:); f(end,end) = f(1,1); 
  else  % interface ? - TODO make seal up poss w/ Bloch phases
    f = fr;
  end
  if fig, figure; end   %else, hold on; end
  h = surf(X,Y,Z, f); set(h,'ambientstrength',0.7,'facelighting','gouraud')
  shading interp; xlabel('x'); ylabel('y'); zlabel('z');
  if fig, lightangle(45,0); end
  axis equal vis3d; 
end
colorbar; colormap(jet(256)); caxis(cs);

%%%%%%%%%%%%
function test_showsurffunc
type = 'torus'; so.a = 1; so.b = 0.5; [s N] = create_panels(type,so);
x = getallnodes(s);
f = sin(1.7*x(1,:)-0.5*x(2,:)+0.9*x(3,:))';    % col vec func
showsurffunc(s, f);
