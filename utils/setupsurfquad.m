function s = setupsurfquad(s,N)
% SETUPSURFQUAD  Set up quadrature nodes on given analytic surf of general type.
%
% s = setupsurfquad(s,N) where N = [Nu Nv], adds spectral surface quadrature
%  nodes to an analytic 3D surface of general type. Types supported:
%   s.type = 't' : torus-like, double-PTR of N=Nu*Nv nodes covering (u,v) in
%                  [0,2pi)^2.
%   s.type = 's' : sphere-like, Nv Gauss-Legendre nodes in v direction, with a
%                  PTR ring in u for each v-node. Variable number of nodes on
%                  each ring is scaled as expected for sphere. Input Nu then
%                  controls overall scaling of numbers on each ring, with
%                  N ~ (2/pi)*Nu*Nv. Note (u,v) is in [0,2pi)x[-1,1].
%
%  Note: is a wrapper to setupdoubleptr or setupspherequad appropriately.
%
%  For definitions of surface quadrature fields, and (u,v) coords of nodes,
%  see: setupdoubleptr, setupspherequad

% Barnett 9/5/19
if s.topo=='t'      % torus type
  s = setupdoubleptr(s,N);
elseif s.topo=='s'  % sphere type
  s = setupspherequad(s,N);
end
