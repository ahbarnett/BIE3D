function s = setupspherequad(s,Ns,o)
% SETUPSPHEREQUAD.  Add PTRxGL node info to an analytic 3D sphere-like surface
%
% s = setupspherequad(s,Ns) where Ns = [Nu Nv] adds a spectral quadrature,
%  periodic trapezoid rule on each ring in u, for each node of a Gauss-Legendre
%  rule in v. The surface must be of type topo='s'. It is defined by functions:
%    s.Z(u,v) - the defining global map from u,v in [0,2pi)x[-1,1] to R3,
%    s.Zu(u,v) and s.Zv(u,v) - its partials.
%  s.Nu becomes the list of numbers of nodes on each ring (even integers), which
%    are bounded by roughly the requested Nu, and are even integers.
%  s.Nv is Nv rounded up to an even integer.
%  If Ns absent or empty, default used.
%
% Note: for unit-aspect (sphere), Nu=2*Nv gives quasi-uniform nodes (h_1~h_2
%  everywhere, apart from rings near poles).
%
% s = setupspherequad(s,Ns,o) also controls opts:
%   o.tensor = 1 (tensor product quadrature, uses pi/2 factor more nodes than
%                 necessary, but simpler), 0 (Nu scaled with elevation, as for a
%                 sphere; more bookkeeping but quasi-optimal).
%            Note: 0 is recommended, and 1 is somewhat obsolete.
%   o.extraunodes, o.minunodes - settings affecting small rings for o.tensor=0.
%
% Output struct s has new node-info fields including:
%  s.x - (3*N) node locs in R3.
%  s.nx - (3*N) unit outward normals.
%  s.xu, s.xv - (3*N) tangential vectors dx/du, dx/dv, for x(u,v) in R3.
%  s.sp - (1*N) "speeds", ie det(Jacobean) from (u,v) -> dS.  sp = cross(xu,xv).
%  s.w - (1*N) weights (including speed).
%  s.N - total number of nodes, typically ~ (2/pi)*Nu*Nv (tensor=0 case),
%        or Nu*Nv (tensor=1).
%  s.Nv - number of v nodes.
%  s.Nu - array of numbers of u-nodes at each v (tensor=0), or number u nodes.
%  s.v - 1d grid of v nodes.
%  s.vw - corresp G-L weights for the nodes s.v for smooth quadr on [-1,1].
%  s.hmin, s.hmax - (1*N) local max, min node spacings.
%
% The node ordering is fast along the incr u direction, slow along incr v.
%
% Is the sphere analog of: setupdoubleptr

% to do: add spectral diff to get s.xu, etc, from only s.x and s.Nu ?

% Barnett 9/3/19
if nargin==0, test_setupspherequad; return; end
if nargin<2 || isempty(Ns), Ns = [60,30]; end
if nargin<3, o=[]; end
if ~isfield(o,'tensor'), o.tensor=0; end    % default
if ~isfield(o,'extraunodes'), o.extraunodes=0; end  % how many above sin(th) val
if ~isfield(o,'minunodes'), o.minunodes=8; end     % handles small rings

Nv=round(Ns(2)/2)*2;    % make even
[s.v s.vw] = gauss(Nv);  % u,v grids (0-offset in u)
if o.tensor
  Nu=round(Ns(1)/2)*2;    % make even
  N = Nu*Nv;
  s.u = (0:Nu-1)'/Nu*2*pi;      % 0-indexed
  [uu vv] = ndgrid(s.u, s.v);   % tensor prod
  uu = uu(:)'; vv = vv(:)';     % since s.Z, cross, etc, row-vectorize only
else
  sintheta = sqrt(1-s.v(:).^2); % for the round sphere only
  Nu = round(o.extraunodes + Ns(1)*sintheta);    % list, scale ~ sin(theta)
  Nu = sqrt(o.minunodes^2 + Nu.^2);     % soft version of max(o.minunodes,Nu)
  Nu = round(Nu/2)*2;     % make all even
  N = sum(Nu);
  uu = nan(1,N); vv=uu;  % alloc output arrays
  offset = 0;
  for i=1:Nv             % loop over const-u (elevation) rings
    Ni = Nu(i);          % this Nu
    uu(offset+(1:Ni)) = (0:Ni-1)'/Ni*2*pi;      % 0-indexed
    vv(offset+(1:Ni)) = s.v(i);
    offset = offset+Nu(i);
  end
end
s.x = s.Z(uu,vv);           % eval all nodes
s.N = N; s.Nv=Nv; s.Nu=Nu;  % write crap out
if isfield(s,'Zu')          % (false allows a nodes-only surface to be formed)
  s.xu = s.Zu(uu,vv);
  s.xv = s.Zv(uu,vv);
  s.nx = cross(s.xu, s.xv);       % outward normal
  s.sp = sqrt(sum(s.nx.^2,1));    % "speeds"
  s.nx = bsxfun(@times,s.nx,1./s.sp);   % unit normal
  if o.tensor
    s.w = (2*pi/Nu)*bsxfun(@times,s.vw, reshape(s.sp,[Nu Nv]));   % quad weights incl "speed"
    s.w = s.w(:)';                  % make row
    hu = sqrt(sum(s.xu.^2,1))*(2*pi/Nu);  % local node sep in u-direc (h_1)
    hv = sqrt(sum(s.xv.^2,1))*(2*pi/Nv);  % local node sep in v-direc (h_2)
  else
    [s.w,hu,hv] = deal(nan(1,N));  % build weights, etc, for each ring in turn
    offset = 0;
    for i=1:Nv
      Ni = Nu(i);          % this Nu
      ii = offset+(1:Ni);  % inds 
      s.w(ii) = (2*pi/Ni) * s.vw(i) * s.sp(ii);
      hu(ii) = sqrt(sum(s.xu(:,ii).^2,1)) * (2*pi/Ni);  % local h_1
      hv(ii) = sqrt(sum(s.xv(:,ii).^2,1)) * s.vw(i);    % local h_2
      offset = offset+Nu(i);
    end    
  end
  s.hmax = max([hu;hv]);        % elementwise, ok when not skew
  s.hmin = min([hu;hv]);        % "
end

%%%%%%%%%
function test_setupspherequad     % tests with sphere, ellipsoid; two quad types
figure;
for shape = 0:1
  if shape==0, disp('sphere'); s = ellipsoid(1,1,1);
  else, disp('ellipsoid'); s = ellipsoid(0.9,1.4,2); end
  for tensor=[1 0]
    o.tensor = tensor; fprintf('tensor=%d:\n',tensor)
    s = setupspherequad(s,[],o);
    tsubplot(2,2,1+shape+2*tensor);
    showsurf(s,'b',struct('alpha',0.2)); lightangle(45,0);
    if shape==0, fprintf('default N=[%d,%d]: surf area err = %.3g\n',max(s.Nu),s.Nv,sum(s.w) - 4*pi); end
    disp('Gauss'' law flux convergence...')
    zo = [0.3; -0.2; 0.5];    % src pt, must be inside the shape
    hold on; plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
    for Nu = 20:20:80, Nv = round(0.5*Nu);       % tie minor discr to major azim
      s = setupspherequad(s,[Nu,Nv],o);
      d = bsxfun(@minus,s.x,zo); r = sqrt(sum(d.^2,1));
      ddotn = sum(d.*s.nx,1);
      flux = sum(s.w.*ddotn./r.^3)/(4*pi);      % surf flux of monopole at zo
      fprintf('N=[%d,%d], N=%d:  \terr = %.3g\n',Nu,Nv,s.N,flux - 1.0);
    end
  end
end
set(gcf,'name','sphere & ellipsoid: quasi-uniform vs tensor-prod quads');
