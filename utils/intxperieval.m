function f = intxperieval(fun,n,y)
% INTXPERIEVAL  Eval func handle (x,y) on periodic-on-rings grid
%
% f = intxperieval(fun,n,y) returns row vector of values of fun(x,y) for nodes
%  (x,y) on periodic 1d-grids lying on constant y=y(i) lines, with n(i)
%  grid points for the ith ring. The domain is [0,2pi)x[-1,1]. The nodes are
%  as described in intxperiinterp.

% split out, 9/5/19
f = nan(1,sum(n));
off = 0;
for i=1:numel(n)
  x = 2*pi*(0:n(i)-1)/n(i);      % note 0-offset!
  f(off+(1:n(i))) = fun(x,y(i));
  off = off + n(i);
end
