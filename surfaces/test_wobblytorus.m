% test a modulated (wobbly) torus. Barnett 10/10/18.
clear
type = 'torus';
so.a = 1.0;   % major radius
b = 0.5;    % mean poloidal radius
wc = 0.1;  % surf modulation ampl
wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
wn = 5;   % # wobbles in toroidal, major
f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
so.b = {f,ft,fp};     % pass in instead of b param
so.p = 6;
o = [];
[pan N] = create_panels(type,so,o);
x = getallnodes(pan);

%figure; plot3(x(1,:),x(2,:),x(3,:),'.','markersize',1); axis equal vis3d
f=0*x(1,:); h=showsurffunc(pan, f);
