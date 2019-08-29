function fff = cruller(b,wc,wa,wb)
% CRULLER - modulation function handles for cruller (twisted wobbly torus)
%
%  fff = cruller(b,wc,wa,wb) returns modulation functions for a cruller, namely
%   cell array of function handles {f(u,v), f_u(u,v), f_v(u,v)} for minor
%   radius-modulation. There are wa (wb) wiggles per major (minor) revolution.
% Inputs:
%  b - baseline minor radius (>0)
%  wc - surf modulation ampl (0 for plain torus; <0 for sin rather than cos)
%  wa - # wobbles around major radius, toroidal
%  wb - # wobbles around minor radius, poloidal
%
% Example: s = setup_torus_doubleptr(1.0,cruller(0.5,0.1,5,3));
%  creates a global discretization of the standard cruller from t-domain BIE
%  paper.
%
% See also: modulatedtorus, setup_torus_doubleptr, create_panels

% Barnett 8/21/19
if wc>=0
  f = @(p,t)  b + wc*cos(wb*t+wa*p);         % wobble func, then its partials...
  fp = @(p,t) -wc*wa*sin(wb*t+wa*p);         % ordering: p=phi=u, then t=theta=v
  ft = @(p,t) -wc*wb*sin(wb*t+wa*p);
else                                         % phase shift by pi/2
  wc = abs(wc);
  f = @(p,t)  b + wc*sin(wb*t+wa*p);         % wobble func, then its partials...
  fp = @(p,t) wc*wa*cos(wb*t+wa*p);          % ordering: p=phi=u, then t=theta=v
  ft = @(p,t) wc*wb*cos(wb*t+wa*p);
end
fff = {f,fp,ft};
