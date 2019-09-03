function s = setupspherequad(s,Ns)
% SETUPSPHEREQUAD.  Add PTRxGL node info to an analytic 3D sphere-like surface
%
% s = setupspherequad(s,Ns) where Ns = [Nu Nv] adds a Nu*Nv quadrature, periodic
%  trapezoid rule on each ring in u, for each node of a Gauss-Legendre rule in
%  v. It is not in general a product quadrature. The surface must be of
%  type topo='s'. It is defined by functions
%  s.Z(u,v) - the defining global map from u,v in [0,2pi)x[-1,1] to R3,
%  s.Zu(u,v) and s.Zv(u,v) - its partials.
%  Nu and Nv are rounded up to even numbers.
%
% Output struct s has new node-info fields including:
%  s.x - (3*N) node locs in R3
%  s.nx - (3*N) unit outward normals
%  s.xu, s.xv - (3*N) tangential vectors dx/du, dx/dv, for x(u,v) in R3.
%  s.sp - (1*N) "speeds", ie det(Jacobean) from (u,v) -> dS.  sp = cross(xu,xv)
%  s.w - (1*N) weights (including speed)
%  s.u, s.v - u and v param grids
%  s.N, S.Nu, s.Nv - numbers of nodes, N = Nu*Nv.
% The node ordering is fast along the incr u direction, slow along incr v.
%
% Is the sphere analog of: setupdoubleptr
%
% Note: product quadr for now (pi/2 not as a efficient at quasi-unif)

% Barnett 9/3/19
if nargin==0, test_setupspherequad; return; end
if nargin<2, Ns = [60,30]; end
Nu=ceil(Ns(1)/2)*2; Nv=ceil(Ns(2)/2)*2;    % make even
s.Nu=Nu; s.Nv=Nv; s.N = Nu*Nv;
s.u = (0:Nu-1)'/Nu*2*pi; [s.v s.vw] = gauss(Nv);  % u,v grids (0-offset in u)
[uu vv] = ndgrid(s.u, s.v); % product for now
uu = uu(:)'; vv = vv(:)';   % since s.Z, cross, etc only vectorize along rows 
s.x = s.Z(uu,vv);
if isfield(s,'Zu')          % allows a nodes-only surface to be formed
  s.xu = s.Zu(uu,vv);
  s.xv = s.Zv(uu,vv);
  s.nx = cross(s.xu, s.xv);       % outward normal
  s.sp = sqrt(sum(s.nx.^2,1));    % "speeds"
  s.nx = bsxfun(@times,s.nx,1./s.sp);   % unit normal
  s.w = (2*pi/Nu)*bsxfun(@times,s.vw, reshape(s.sp,[Nu Nv]));   % quad weights incl speed
  s.w = s.w(:)';                  % make row
end

%%%%%%%%%
function test_setupspherequad     % tests with ellipsoid
for shape = 0:1
  if shape==0, disp('sphere'); s = ellipsoid(1,1);
  else, disp('ellipsoid'); s = ellipsoid(1.5,2); end
  s = setupspherequad(s);
  figure; showsurf(s,'b',struct('alpha',0.2)); lightangle(45,0);
  if shape==0, fprintf('default N=[%d,%d]: surf area err = %.3g\n',s.Nu,s.Nv,sum(s.w) - 4*pi); end
  disp('Gauss'' law flux convergence...')
  zo = [0.5; -0.2; 0.1];    % src pt, must be inside the shape
  hold on; plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
  for Nu = 20:20:80, Nv = ceil(0.5*Nu);       % tie minor discr to major azim
    s = setupspherequad(s,[Nu,Nv]);
    d = bsxfun(@minus,s.x,zo); r = sqrt(sum(d.^2,1));
    ddotn = sum(d.*s.nx,1);
    flux = sum(s.w.*ddotn./r.^3)/(4*pi);      % surf flux of monopole at zo
    fprintf('N=[%d,%d]:  \terr = %.3g\n',Nu,Nv,flux - 1.0);
  end
end
