function s = modulatedtorus(a,b)
% MODULATEDTORUS.  Analytic global description of radius-modulated torus surface
%
% s = modulatedtorus(a,b) returns a surface struct with function handles
%  Z, Zu, Zv : [0,2pi)^2 -> R^3,
%  describing analytically a torus in 3D, with major radius a, minor radius b.
%  Z is the map, and Zu, Zv its partials.
%  1st coord (u) goes CCW (viewed from above) around big circle, toroidal
%  2nd coord (v) takes up around little circle, poloidal.
%  Zu cross Zv points outwards.
%  The function handles vectorize over rows (only) of input values.
%
% If b is a cell array of 3 function handles, they will be interpreted
%  as the (minor) radial function {f(u,v),f_u(u,v),f_v(u,v)}, where the
%  user-supplied partials must be correct.
%
% Based on the obsolete torusparam.m, but correcting the flipped f(.,.) order

% Barnett 8/21/19
if ~iscell(b)        % plain torus.  p=phi=u, t=theta=v
  s.Z = @(p,t) [(a+b*cos(t)).*cos(p); (a+b*cos(t)).*sin(p); b*sin(t)];
  s.Zu = @(p,t) [-(a+b*cos(t)).*sin(p); (a+b*cos(t)).*cos(p); 0*p];
  s.Zv = @(p,t) [-b*sin(t).*cos(p); -b*sin(t).*sin(p); b*cos(t)];
else                 % modulated torus
  f = b{1}; fp = b{2}; ft = b{3};
  s.Z = @(p,t) [(a+f(p,t).*cos(t)).*cos(p); (a+f(p,t).*cos(t)).*sin(p); f(p,t).*sin(t)];
  s.Zu = @(p,t) [-(a+f(p,t).*cos(t)).*sin(p)+fp(p,t).*cos(t).*cos(p); (a+f(p,t).*cos(t)).*cos(p) + fp(p,t).*cos(t).*sin(p); fp(p,t).*sin(t)];
  s.Zv = @(p,t) [(-f(p,t).*sin(t) + ft(p,t).*cos(t)).*cos(p); (-f(p,t).*sin(t)+ft(p,t).*cos(t)).*sin(p); f(p,t).*cos(t)+ft(p,t).*sin(t)];
end
