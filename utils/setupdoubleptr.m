function s = setupdoubleptr(s,Ns)
% SETUPDOUBLEPTR.  Add double-PTR node info to an analytic 3D torus-like surface
%
% s = setupdoubleptr(s,N) where N = [Nu Nv] adds a Nu*Nv double periodic
%  trapezoid rule quadrature to a torus-like global surface defined by funcions
%  s.Z(u,v) - the defining global map from u,v in [0,2pi) to R3,
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
% The node ordering is fast along the u direction, slow along v.
%
%  If the partials are not present, nodes s.x only are added, for now
%  *** to do: add double spectral diff here.
%
% It is the 3D analog of setupquad in BIE2D.

% Barnett 8/21/19
if nargin==0, test_setupdoubleptr; return; end
if nargin<2, Ns = [60,30]; end
Nu=ceil(Ns(1)/2)*2; Nv=ceil(Ns(2)/2)*2;    % make even
s.Nu=Nu; s.Nv=Nv; s.N = Nu*Nv;
s.u = (0:Nu-1)'/Nu*2*pi; s.v = (0:Nv-1)'/Nv*2*pi;  % u,v grids, 0-offset
[uu vv] = ndgrid(s.u, s.v);
uu = uu(:)'; vv = vv(:)';   % since s.Z, cross, etc only vectorize along rows 
s.x = s.Z(uu,vv);
if isfield(s,'Zu')          % allows a nodes-only surface to be formed
  s.xu = s.Zu(uu,vv);
  s.xv = s.Zv(uu,vv);
  s.nx = cross(s.xu, s.xv);       % outward normal
  s.sp = sqrt(sum(s.nx.^2,1));    % "speeds"
  s.nx = bsxfun(@times,s.nx,1./s.sp);   % unit normal
  s.w = (2*pi/Nu)*(2*pi/Nv) * s.sp;   % quad weights incl speed
end
  
%%%%%%%%
function test_setupdoubleptr
a=1.0; b= 0.5; s = modulatedtorus(a,b);
s = setupdoubleptr(s);
o=[]; o.alpha = 0.2; figure; showsurf(s,'b',o); lightangle(45,0);
