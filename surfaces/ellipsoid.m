function s = ellipsoid(a,b)
% ELLIPSOID.  Analytic global ellipsoid aligned with axes
%
% s = ellipsoid(a,b) makes surface struct with only analytic function handles
%  Z, Zu, Zv : [0,2pi)x[-1,1] -> R^3, describing a coord-aligned ellipsoid
%  with semiaxes (1,a,b).  u is phi (azimuthal) and v is z (elevation).
%  Vectorizes over rows for each of u,v.
%
% For test run: setupspherequad

% Barnett 9/3/19
s.Z = @(u,v) [cos(u).*sqrt(1-v.^2); a*sin(u).*sqrt(1-v.^2); b*v];
s.Zu = @(u,v) [-sin(u).*sqrt(1-v.^2); a*cos(u).*sqrt(1-v.^2); zeros(size(v))];
s.Zv = @(u,v) [cos(u).*(-v)./sqrt(1-v.^2); a*sin(u).*(-v)./sqrt(1-v.^2); b*ones(size(v))];
s.topo = 's';   % "sphere" type
