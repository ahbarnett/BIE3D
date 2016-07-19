function [x nx sp dp dt] = torusparam(a,b,p,t,o)
% TORUSPARAM - return coords, unit normal and speed at params (phi,theta)
%
% [x nx sp dp dt] = torusparam(a,b,p,t)  p,t should be row vectors
%  1st coord (phi) goes CCW (viewed from above) around big circle
%  2nd coord (theta) takes up around little circle
% (dp,dt,normal) form a RH coord sys.
%
% sp = dp * dt
% dt is magnitude of dr/dt vector; dp mag of dr/dp vector.
%
% [x nx sp dp dt] = torusparam(a,b,p,t,o) allows opts eg
% o.f, o.fp, o.ft handles to modulation func (poloidal radius), and partials

% adapted from Barnett 2013 except made outward normal dphi cross dtheta

if nargin>4       % wobbly torus
  f = o.f; ft = o.ft; fp = o.fp;
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
