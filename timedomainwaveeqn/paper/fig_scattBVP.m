% t-domain wave eqn BIE: scattering (exterior Dirichlet BVP).
% paper fig version, Barnett 1/9/19, based on fig_toggleab.m
% fsparse crashes for big probs (on desktop not laptop)
%rmpath ~/matlab/stenglib/Fast
clear
dt =   0.1;   % timestep
Ttot = 15.0;     % total time to evolve
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
so.np = 6; so.mp=so.np*2/3;    % panel numbers in major and minor direction
o.nr = 2*o.p; o.nt = 2*o.nr;  % aux quad orders: radial and angular nodes
% Dirichlet data for BVP: t-dep pulse func
t0=0; s0=0.5; T = @(t) 50*exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
xs = [4.5;0.1;-0.2];   % src pt for data, now outside

[s N] = create_panels('torus',so,o);     % surf
[x nx w] = getallnodes(s);
gdata = @(t) data_ptsrc(xs,T,Tt,t,x);    % Dirichlet data *func*
wpred = extrap(m);          % extrapolation row vector
fprintf('surf quad order p=%d, %dx%d panels, N=%d\n',o.p,so.np,so.mp,N)

distmax = 6.0;                   % largest dist from anything to anything
n = ceil(distmax/dt);
s = add_panels_auxquad(s,o);
Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
t0=tic;
[Starg,Dtarg,Sdottarg] = tdSDinterpmats_panels(s,s,Linfo,struct('n',n,'dt',dt,'m',m));
fprintf('tot S,D,Sdot mat build (each %dx%d, nnz=%d): %.3g s\n',size(Starg,1),size(Starg,2),nnz(Starg),toc(t0))

% GRF test: eval sig,tau on {n history grid} x {N bdry nodes}; time fast, N slow
ttarg=6.0; tt = ttarg + dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,T,Tt,ttt,xx,nxx);       % (output ft unused)
sighist = -fn; tauhist = f;  % col vecs, ext wave eqn GRF: u = D.u - S.un
%f = reshape(f,[n N]); figure;
%oo.nofig=1; oo.sc=1; for i=1:n, showsurffunc(s,f(i,:),oo); drawnow, end

% set up vectors which compute potential u(t,x) at >1 targs x, from dens hist...
t.N = 100; t.x = [1.0;0.2;0.7] + [0;0;3]*(1:t.N)/t.N;  % targs
showsurffunc(s,0*x(1,:)); hold on; plot3(xs(1),xs(2),xs(3),'r.');
plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.');

tret = -dists(t.x,x);  % retarded times of surf nodes rel to ttarg, trg loc fast
[jmax,jmin,a,ap] = interpmat(tret,dt,m);   % Tom's coeffs, 1 row per tret entry
joff = jmin+n-1;         % padding on the ancient side
if joff<0, error('u_test eval vec requesting too ancient history!'); end
a = [zeros(N*t.N,joff), a, zeros(N*t.N,-jmax)];  % pad to width n
ap = [zeros(N*t.N,joff), ap, zeros(N*t.N,-jmax)];
[S D Dp] = tdSDmats(t.x,x,nx,w);  % each is NtxN (Nt = N.t = # targ)
Stest = a.*repmat(S(:),[1 n]);  % coeff vectors packed as (Nt*N) x n mats
Dtest = a.*repmat(D(:),[1 n]) + ap.*repmat(Dp(:),[1 n]);
Sdottest = ap.*repmat(S(:),[1 n]);
Stest = reshape(permute(reshape(Stest,[t.N N n]),[1 3 2]),[t.N N*n]);
Dtest = reshape(permute(reshape(Dtest,[t.N N n]),[1 3 2]),[t.N N*n]);
Sdottest = reshape(permute(reshape(Sdottest,[t.N N n]),[1 3 2]),[t.N N*n]);
utest = Stest*sighist + Dtest*tauhist;   % two terms in eval GRF
fprintf('GRF null-field test, max at all %d targs: %.3g\n', t.N, max(abs(utest)))

al=1; be=2;   % Representation is u = D.mu + be.S.mu + al.S.mudot:
Rtarg = Dtarg + al*Sdottarg + be*Starg;  % for history application of ret BIEs
Rtest = Dtest + al*Sdottest + be*Stest;  % for the test pt eval
mo = []; mo.verb = 1; mo.shift=corrshift;  % opts for timestepping
[tj u rhsnrm gnrm munrm muall] = tmarch(dt,Ttot,predcorr,gdata,Rtarg,Rtest,wpred,mo);  % do time steps

figure; imagesc(tj,t.x(3,:),u); colorbar;
ylabel('z coord of targ'); xlabel('t time');



