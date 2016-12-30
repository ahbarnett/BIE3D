function [pan N] = create_panels(shape,so,o)
% CREATE_PANELS  Build panels for various 3D surface geometries
%
% [pan N] = create_panels(shape,shapeopts,convopts)
%  Creates struct array of panels with their quadrature nodes for unknowns
%
% eg shape = 'torus':    shapeopts has fields a (major radius), b (minor rad),
%                        mp (# panels small circ), np (# pan big circ).
%                        mp * np panels are ordered in incr theta (i) fast,
%                        incr phi (j) slow  (ie flipud of matrix ordering).
%                        Dof ordering within panels is same (y fast, x slow)
%
% Output
%  N : total number of unknowns
%  pan : cell array with each element having fields
%        (n = p^2 is # nodes per panel):
%        inds - indices from full N list belonging in this panel
%        t - (2 by n) chart parameters of smooth nodes
%        x - (3 by n) 3D coords of smooth nodes
%        nx - (3 by n) 3D coords of normal derivs at smooth nodes
%        w - (1 by n) smooth quad weights including speed factors
%        [sp - (1 by n) speed factors at nodes (Jacobean area factor)]
%        chart - parametrization function for panels, specifically:
%                input 2-by-N pts in [-1,1]^2, output [x n_x speed] where
%                x and nx are 3-by-N (pts and normals in R3), speed is 1-by-N.
%                (note speed is relative to global param of shape)
%        spfac - (scalar) size factor converting to speed in chart coords.
%        nei - indices in panel list of (usually 8) neighbor panels, in a
%              standard ordering (must match that in set_auxinterp).
%        neix1 - (3 by 8) 3D coords of first node in each neighbor panel,
%              in std ordering, useful to identify neighbors w/o global index
%              (used by: relatedpanel.m)
%  In addition pan{1} contains the following fields needed to specify geometry:
%        mp, np - numbers of vertical and horizontal panels.
%        topo - topology, eg 'torus'

% Barnett 7/11/15. restarted 7/14/16, swapped np & mp 7/18/16
% 12/29/30 neix1.

% future ideas:
% *** why not make just the surf pts be input, and derive all geom
% from that? But, let the q used for that differ from the p used for
% density rep.
%   *** Use bkwd stable monomial solve to get the local charts for
%  r, r_u, r_v, n, etc ?

if nargin==0, test_create_panels; return; end
if nargin<3, o=[]; end
if isfield(o,'p'), p = o.p; else p = 8; end  % # nodes per side of panel (order)

if strcmp(shape,'torus') 
  if isfield(so,'mp'), mp = so.mp; else mp = 8; end  % # pans on small circle
  if isfield(so,'np'), np = so.np; else np = 12; end  % # pans on big circle
  tpansiz1 = pi/np; tpansiz2 = pi/mp;     % panel half-sizes in parameter
  npan = mp*np; N = npan*p^2;       % same order all panels, total # nodes
  panind = reshape(1:npan,[mp np]); % 2d grid of panel numbers (sets dof order)
  pan = cell(npan,1);
  for j=1:np, for i=1:mp, k=panind(i,j);     % loop over panels (in any order)
      pan{k}.psiz = [tpansiz1,tpansiz2];     % param sizes
      pan{k}.inds = p^2*(k-1)+(1:p^2);     % set up indices of unknowns
      pan{k}.chart = @(t) torusparam(so.a,so.b,(t(1,:)+2*j)*tpansiz1,(t(2,:)+2*i)*tpansiz2);   % note speed is relative to original params; note i up, j across
      pan{k}.spfac = tpansiz1*tpansiz2;  % factor to convert speed to chart z
      nei = panind(mod([i-2 i-1 i],mp)+1,mod([j-2 j-1 j],np)+1); % ordering key
      nei = nei(nei~=k);    % remove self from list, to leave 8
      pan{k}.nei = nei(:);
    end,end
    pan{1}.mp=mp; pan{1}.np=np; pan{1}.topo = 'torus'; pan{1}.p = p;
end
pan = panel_smooth_quad(pan,p);  % does all panels
% add neix1 tags...
for k=1:numel(pan), nei=pan{k}.nei;
  x = []; for n=1:numel(nei), x = [x pan{nei(n)}.x(:,1)]; end
  pan{k}.neix1 = x;
end
%%%

function pan = panel_smooth_quad(pan,p)   % setup G-L tensor native quadr
[x w] = gauss(p);   % Gauss-Legendre nodes & weights on [-1,1]
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
ww = w(:)*w; ww = ww(:)';                     % 1*p^2 weights
for m=1:numel(pan)
  [pan{m}.x pan{m}.nx sp] = pan{m}.chart(xx);  % get 3 outputs
  pan{m}.w = sp.*ww*pan{m}.spfac;
  pan{m}.N = p*p;
  pan{m}.t = xx;          % same for all panels, for now
end

%%%%%%%%%%%%%%%%%%%%%%
function test_create_panels  % self-test
type = 'torus', so.a = 1; so.b = 0.5; o = [];
[pan N] = create_panels(type,so,o);

% plot it...
figure; x = getallnodes(pan);
plot3(x(1,:),x(2,:),x(3,:),'.','markersize',1); axis equal vis3d
k = 57;    % random panel (conveniently faces viewer for np=8, mp=12)
hold on; for i=1:8, knei=pan{k}.nei(i);
  x = pan{knei}.x; plot3(x(1,:),x(2,:),x(3,:),'g.','markersize',5);
  l = mean(x,2); text(l(1),l(2),l(3),sprintf('%d(%d)',knei,i+1),'fontsize',20)
end
hold on; x = pan{k}.x; plot3(x(1,:),x(2,:),x(3,:),'r.','markersize',10); % self
l = mean(x,2); text(l(1),l(2),l(3),sprintf('%d(%d)',k,1),'fontsize',20,'color',[1 0 0])
l = 0.1; nx = pan{k}.nx; plot3([x(1,:);x(1,:)+l*nx(1,:)],[x(2,:);x(2,:)+l*nx(2,:)],[x(3,:);x(3,:)+l*nx(3,:)],'k-');  % normals
text(x(1,1),x(2,1),x(3,1),'dof 1'); text(x(1,2),x(2,2),x(3,2),'dof 2');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
pan{k}   % look inside struct
%pan{k}.nei

disp('surface area test:')
np = 3; so.np = np; so.mp = ceil(1.5*np); 
[pan N] = create_panels(type,so,o);
[~,~,w] = getallnodes(pan);    % get all weights
fprintf('np=%d: area err = %.3g\n',np,sum(w) - 2*pi^2)

disp('Gauss law flux test:')
zo = [0.9; -0.2; 0.1];    % src pt, must be inside the shape
plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
for mp = 4:2:12, so.mp = mp; so.np = ceil(1.5*mp); 
  [pan N] = create_panels(type,so,o);
  [x,nx,w] = getallnodes(pan);
  d = bsxfun(@minus,x,zo); r = sqrt(sum(d.^2,1));
  ddotn = sum(d.*nx,1);
  fprintf('mp=%d (N=%d): flux err = %.3g\n',mp,N,sum(w.*ddotn./r.^3)/(4*pi) - 1.0)
end
