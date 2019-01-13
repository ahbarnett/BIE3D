% Show sparsity patterns of A^r_{ij}. Barnett 12/21/18.

clear
%%%%  all params...
dt = 0.05;   % timestep
Ttot = 27.0;     % total time to evolve
m = 4;      % control time interp order (order actually m+2)
predcorr = 8;   % <0 for impl;  0,1,2, for pred-corr with that many corr steps
corrshift = 0.25;    % only helps for large dt
wobbly = 0;   % 0: torus, 1: cruller.
so.a=1; b=0.5; o.p = 6;  % torus shape (a,b);  p = G-L nodes per panel side.
if ~wobbly  % plain torus
  so.b =b;
else        % cruller
  wc = 0.1;  % surf modulation ampl
  wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
  wn = 5;   % # wobbles in toroidal, major
  f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
  ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
  so.b = {f,ft,fp};     % pass in instead of b param
end
so.np = 12; so.mp=so.np*2/3;    % panel numbers in major and minor direction.
o.nr = 2*o.p; o.nt = 2*o.nr;  % aux quad orders: radial and angular nodes
% Dirichlet data for BVP: t-dep pulse func for interior pt src... (max val ~1)
t0=6; s0=1.0; T = @(t) 5*exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
% (t0/s0 = 6 gives 1e-8 of start-up error if no time delay from src to surf)
%eps = 1e-5; tt = 4.3; fprintf('check Tt vs T: %.3g\n',(T(tt+eps)-T(tt-eps))/(2*eps) - Tt(tt)), clear tt
%%%%

[s N] = create_panels('torus',so,o);     % surf
[x nx w] = getallnodes(s);
fprintf('surf quad order p=%d, %dx%d panels, N=%d\n',o.p,so.np,so.mp,N)
distmax = 6.0;                   % largest dist from anything to anything
n = ceil(distmax/dt);
s = add_panels_auxquad(s,o);
Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
t0=tic;
[Starg,Dtarg,Sdottarg] = tdSDinterpmats_panels(s,s,Linfo,struct('n',n,'dt',dt,'m',m));
fprintf('tot S,D,Sdot mat build (each %dx%d, nnz=%d): %.3g s\n',size(Starg,1),size(Starg,2),nnz(Starg),toc(t0))

R = Starg+Dtarg; % dummy
figure; rr=0:5:60; nr=numel(rr);     % plot them .............
for s=1:nr, r = rr(s); 
  h=axes('position',[(0.05+s-1)/nr, 0.05, 0.9/nr, 0.9]);
  spy(R(:,n-r:n:end)); axis off; title(sprintf('r=%d',r));
end
set(gcf,'paperposition',[0 0 10 1.5]);
print -depsc2 Arij.eps

N
'sparsity =', nnz(R)/numel(R)

