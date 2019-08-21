function b = cruller(b,wc,wa,wb)
% CRULLER - modify torus b param to make it a cruller (twisted wobbly torus)
%
%  b = cruller(b,wc,wm,wn) returns 2nd parameter for a torus which is a struct
%   of function handles giving torus radius-modulation, for the generalized
%   torus.
%  b - baseline minor radius
%  wc - surf modulation ampl (0 for plain torus)
%  wa - # wobbles around major radius, toroidal
%  wb - # wobbles around minor, poloidal
%
% Example: s = setup_torus_doubleptr(1.0,cruller(0.5,0.1,5,3));
%  creates a global discretization of the standard cruller from t-domain BIE
%  paper.
%
% See also: torusparam, setup_torus_doubleptr, create_panels

% Barnett 8/21/19
f = @(t,p) b + wc*cos(wb*t+wa*p);         % wobble func, then its partials...
ft = @(t,p) -wc*wb*sin(wb*t+wa*p); fp = @(t,p) -wc*wa*sin(wb*t+wa*p);
b = {f,ft,fp};     % replaces b param
