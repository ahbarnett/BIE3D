function f = fun2dquadeval(fun,s)
% FUN2DQUADEVAL  Evaluate function handle over param nodes of spectral 2d quadr
%
% f = fun2dquadeval(fun,s) returns a row-vec of values of fun(u,v) for each of
%  (u,v) in the set of nodes in the parameter space for the surface s.
%  s.topo determines the parameter space: 't' torus is [0,2pi)^2, whereas
%  's' sphere is [0,2pi)x[-1,1]. The grid parameters are taken from the surface
%  struct s (see setupdoubleptr or setupspherequad). fun must vectorize along
%  rows; ie, u and v can be row vectors of equal lengths.

% Barnett 9/5/19
if s.topo=='t'      % torus-like
  [uu vv] = ndgrid(s.u, s.v);   % tensor prod
  uu = uu(:)'; vv = vv(:)';
  f = fun(uu,vv);
elseif s.topo=='s'  % sphere-like
  f = intxperieval(fun,s.Nu,s.v);
end
