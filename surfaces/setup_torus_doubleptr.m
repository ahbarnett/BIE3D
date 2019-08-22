function s = setup_torus_doubleptr(a,b,Ns)
% SETUP_TORUS_DOUBLEPTR.  analytic torus-like 3D surf w/ global double-PTR quad
%
% s = setup_torus_doubleptr(a,b,Ns) creates double-PTR global quadrature
%  of a (possibly radius-modulated) torus. See modulatedtorus for geometry (a,b)
%  and setupdoubleptr for Ns. Now only a simple wrapper and test code.

% Simplified from create_panels.m (2015) and torusquad.m (2013).
% Barnett 8/16/19, gutted to be a wrapper, 8/21/19.

if nargin==0, test_setup_torus_doubleptr; return; end
if nargin<3, Ns = [60,30]; end
s = modulatedtorus(a,b);
s = setupdoubleptr(s,Ns);

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
  figure; showsurf(s,'b',struct('alpha',0.2)); lightangle(45,0);
  if shape==0, fprintf('default N=[%d,%d]: surf area err = %.3g\n',s.Nu,s.Nv,sum(s.w) - 2*pi^2), end
  
  disp('Gauss'' law flux convergence...')
  zo = [0.9; -0.2; 0.1];    % src pt, must be inside the shape
  hold on; plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
  for Na = 20:20:80, Nb = ceil(0.5*Na);       % tie minor discr to major
    s = setup_torus_doubleptr(a,b,[Na,Nb]);   % Na aka Nu; Nb aka Nv
    d = bsxfun(@minus,s.x,zo); r = sqrt(sum(d.^2,1));
    ddotn = sum(d.*s.nx,1);
    flux = sum(s.w.*ddotn./r.^3)/(4*pi);      % surf flux of monopole at zo
    fprintf('N=[%d,%d]:  \terr = %.3g\n',Na,Nb,flux - 1.0);
  end
end
