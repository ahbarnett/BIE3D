function [x nx sp] = torusparam(a,b,p,t,o)
% TORUSPARAM - return coords, unit normal and speed at params (phi,theta)
%
% Warning: this is somewhat obsolete, kept only for panels & timedomainwaveeqn.
% Prefer: modulatedtorus (which contains the analytic functions), and
%         setupdoubleptr (which shows how to set up nodes from analytic funcs).
%
% [x nx sp] = torusparam(a,b,p,t).
%  Inputs:
%  p,t are phi and theta, and should be row vectors if multiple points wanted
%  a is major radius
%  b is minor radius.
%    If b is a cell array of 3 function handles, is interpreted as b(t,p),
%    a minor (poloidal) radius function over theta and phi each in [0,2pi),
%    and b_t(t,p), b_p(t,p) its correct partials. Or, see below option.
%
% Geometry:
%  1st coord (phi) goes CCW (viewed from above) around big circle, toroidal
%  2nd coord (theta) takes up around little circle, poloidal
% (dp,dt,outwardsnormal) form a RH coord sys.
%
% [x nx sp] = torusparam(a,b,p,t,o) allows opts eg
% o.f, o.fp, o.ft handles to modulation func (poloidal radius vs theta,phi),
% and its partials, which must be correct.

% adapted from Barnett 2013 except made outward normal dphi cross dtheta.
% expanded modulation interface 10/10/18.

if nargin>4 || iscell(b)       % wobbly torus
  if iscell(b)    % allow a new interface to wobbly modulation
    f = b{1}; ft = b{2}; fp = b{3};
  else
    f = o.f; ft = o.ft; fp = o.fp;
  end
  r = @(t,p) [(a+f(t,p).*cos(t)).*cos(p); (a+f(t,p).*cos(t)).*sin(p); f(t,p).*sin(t)]; % t,p rows
  rt = @(t,p) [(-f(t,p).*sin(t) + ft(t,p).*cos(t)).*cos(p); (-f(t,p).*sin(t)+ft(t,p).*cos(t)).*sin(p); f(t,p).*cos(t)+ft(t,p).*sin(t)]; % partial
  rp = @(t,p) [-(a+f(t,p).*cos(t)).*sin(p)+fp(t,p).*cos(t).*cos(p); (a+f(t,p).*cos(t)).*cos(p) + fp(t,p).*cos(t).*sin(p); fp(t,p).*sin(t)];
else % param...
  r = @(t,p) [(a+b*cos(t)).*cos(p); (a+b*cos(t)).*sin(p); b*sin(t)]; % t,p rows
  rt = @(t,p) [-b*sin(t).*cos(p); -b*sin(t).*sin(p); b*cos(t)]; % partial
  rp = @(t,p) [-(a+b*cos(t)).*sin(p); (a+b*cos(t)).*cos(p); 0*p];
end

x = r(t,p);
rts = rt(t,p); rps = rp(t,p);
dt = sqrt(sum(rts.^2,1));
dp = sqrt(sum(rps.^2,1));    % could reuse this data
nx =  cross(rps, rts);       % outward normal
sp = sqrt(sum(nx.^2,1)); % speeds
nx = nx.*repmat(1./sp, [3 1]);
