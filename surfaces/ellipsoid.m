function s = ellipsoid(a,b,c)
% ELLIPSOID.  Analytic global ellipsoid aligned with axes
%
% s = ellipsoid(a,b,c) makes surface struct with only analytic function handles
%  Z, Zu, Zv : [0,2pi)x[-1,1] -> R^3, describing a coord-aligned ellipsoid
%  with semiaxes (a,b,c).  u is phi (azimuthal) and v is z (elevation).
%  Vectorizes over rows for each of u,v.
%
% For test run: setupspherequad

% Barnett 9/3/19
s.Z = @(u,v) [a*cos(u).*sqrt(1-v.^2); b*sin(u).*sqrt(1-v.^2); c*v];
s.Zu = @(u,v) [a*-sin(u).*sqrt(1-v.^2); b*cos(u).*sqrt(1-v.^2); zeros(size(v))];
s.Zv = @(u,v) [a*cos(u).*(-v)./sqrt(1-v.^2); b*sin(u).*(-v)./sqrt(1-v.^2); c*ones(size(v))];
s.topo = 's';   % "sphere" type
