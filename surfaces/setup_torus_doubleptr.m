function s = setup_torus_doubleptr(a,b,Ns,o)
% SETUP_TORUS_DOUBLEPTR.  analytic torus-like 3D surf, double global PTR quad
%
% s = setup_torus_doubleptr(a,b,Ns,o) creates double-PTR global quadrature
%  of a (possibly radius-modulated) torus.
%
% Definition of coord sys (coming from torusparam.m):
%  1st coord (phi) goes CCW (viewed from above) around big circle, toroidal.
%  2nd coord (theta) takes up around little circle, poloidal.
%  Thus (dp,dt,outwardsnormal) form a RH coord sys.
%
% Inputs:
%  a - major radius
%  b - minor radius, or, a cell array of 3 function handles, interpreted as
%                   { b(theta,phi), b_theta(theta,phi), and b_phi(theta,phi },
%                   a minor radius function over theta and phi in [0,2pi), and
%                   its correct partials. This allows height-modulated torus.
%                   Note that the functions for historical reasons have (t,p)
%                   not the usual (p,t) ordering.
%  Ns - 2-element vector [Na,Nb] where Na=# major nodes, Nb=# minor nodes.
%       They are rounded up to even numbers.
% Output:
%  s global surface struct contains:
%       N (total # nodes), x (3*N) node locs, nx (3*N) unit outwards normals,
%       w (1*N) weights (including speed), etc
%       (This serves as a definition of the global 3D surf struct object.)
%
%       Here for double-PTR, N=Na*Nb, and s.p,s.t (phi,theta 1d grids) output.
%       Nodes are ordered fast in the a (major, phi) dir, going rightwards,
%       slow in the b (minor, theta) dir, going upwards, viewed from (+inf,0,0)
%       ie, "CCW 90deg rot of fortran matrix ordering". Flipped relative to
%       create_panels('torus',...), but more logical.
%
% See also: torusparam.m
%
% Note: this is the 3D analog of BIE2D/utils/setupquad.m (global PTR)

% Simplified from create_panels.m (2015) and torusquad.m (2013)
% Barnett 8/16/19

if nargin==0, test_setup_torus_doubleptr; return; end
if nargin<3, Ns = [60,30]; end
Na=ceil(Ns(1)/2)*2; Nb=ceil(Ns(2)/2)*2;    % make even
s.Na=Na; s.Nb=Nb; s.N = Na*Nb;

s.p = (0:Na-1)'/Na*2*pi; s.t = (0:Nb-1)'/Nb*2*pi;  % phi,theta grids, 0-offset
[pp tt] = ndgrid(s.p, s.t);
[s.x s.nx s.w] = torusparam(a,b,pp(:)',tt(:)');   % speeds go into s.w
s.w = (2*pi/Na)*(2*pi/Nb) * s.w;   % quad weights incl speed


%%%%%%%
function test_setup_torus_doubleptr
for shape = 0:1
  a = 1.0; b = 0.5;       % baseline torus params
  if shape==0
    disp('plain torus double PTR quadr test:')
  else
    disp('cruller double PTR quadr test:')
    b = cruller(b,0.1,5,3);    % replaces b
  end
  s = setup_torus_doubleptr(a,b);
  figure; plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.','markersize',1); axis equal vis3d
  if shape==0, fprintf('default N=[%d,%d]: surf area err = %.3g\n',s.Na,s.Nb,sum(s.w) - 2*pi^2), end
  
  disp('Gauss'' law flux convergence...')
  zo = [0.9; -0.2; 0.1];    % src pt, must be inside the shape
  hold on; plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
  for Na = 20:20:80, Nb = ceil(0.5*Na);      % tie minor discr to major
    s = setup_torus_doubleptr(a,b,[Na,Nb]);
    d = bsxfun(@minus,s.x,zo); r = sqrt(sum(d.^2,1));
    ddotn = sum(d.*s.nx,1);
    flux = sum(s.w.*ddotn./r.^3)/(4*pi);   % surf flux of monopole at zo
    fprintf('N=[%d,%d]:  \terr = %.3g\n',Na,Nb,flux - 1.0);
  end
end
