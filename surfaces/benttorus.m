function fff = benttorus(b,wc,wa)
% BENTTORUS - modulation function handles for bent torus
%
%  fff = benttorus(b,wc,wa) returns modulation functions for a bent torus,
%   cell array of function handles {f(u,v), f_u(u,v), f_v(u,v)} for minor
%   radius-modulation. There are wa wiggles per major revolution.
% Inputs:
%  b - baseline minor radius (>0)
%  wc - surf modulation ampl (0 for plain torus)
%  wa - # wobbles around major radius, toroidal
%
% See also: cruller

% Barnett 8/29/19
f = @(p,t)  b + wc*cos(wa*p).*sin(t);         % wobble func, then its partials...
fp = @(p,t) -wc*wa*sin(wa*p).*sin(t);         % ordering: p=phi=u, then t=theta=v
ft = @(p,t) wc*cos(wa*p).*cos(t);
fff = {f,fp,ft};
