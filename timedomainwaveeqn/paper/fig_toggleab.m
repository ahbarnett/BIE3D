% t-domain wave eqn BIE: known exterior Dirichlet BVP.
% Four different reps, output eps of max norms vs j. Store all results.
% paper fig version, Barnett 12/18/18, based on demo_waveeqn_tmarch_*.m

% fsparse crashes for big probs (on desktop not laptop)
%rmpath ~/matlab/stenglib/Fast
clear

%%%%  all params...
dt = 0.1;   % timestep
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
so.np = 9; so.mp=so.np*2/3;    % panel numbers in major and minor direction.
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

% surf data for GRF test of SDtest vectors only...
w0 = 2.0; TT = @(t) cos(w0*t); TTt = @(t) -w0*sin(w0*t); % data src t-func, tested (don't confuse with emitter func T, below)
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside
% eval sig, tau on {n history grid} x {N bdry nodes}; recall time fast, N slow
tt = dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,TT,TTt,ttt,xx,nxx);       % output ft unused
sighist = -fn; tauhist = f;  % col vecs, ext wave eqn GRF: u = D.u - S.un

% set up vectors which compute potential at fixed u(t,x_ext) from dens hist...
t.N = 1; t.x = [1.3;0.1;0.8];    % single test targ pt, exterior
tret = -dists(t.x,x);     % retarded times of surf nodes rel to test pt
[jmax,jmin,a,ap] = interpmat(tret,dt,m);    % Tom's coeffs (1 row per tret)
joff = jmin+n-1;         % padding on the ancient side
if joff<0, error('u_test eval vec requesting too ancient history!'); end
a = [zeros(N,joff), a, zeros(N,-jmax)];  % pad to width n, preserve sparsity
ap = [zeros(N,joff), ap, zeros(N,-jmax)];
[S D Dp] = tdSDmats(t.x,x,nx,w);  % each is 1xN
Stest = a'.*repmat(S,[n 1]);      % coeff vectors packed as nxN matrices
Dtest = a'.*repmat(D,[n 1]) + ap'.*repmat(Dp,[n 1]);
Sdottest = ap'.*repmat(S,[n 1]);
Stest = Stest(:)'; Dtest = Dtest(:)'; Sdottest = Sdottest(:)'; % as row vecs

utest = Stest*sighist + Dtest*tauhist;   % two terms in eval GRF
uex = data_ptsrc(xs,TT,TTt,0.0,t.x);       % vs exact u at test pt
fprintf('GRF test that u eval vectors (SDtest) work: %.3g\n', utest-uex)
wpred = extrap(m);          % extrapolation row vector
gdata = @(t) data_ptsrc(xs,T,Tt,t,x,nx);    % Dirichlet data func from xs src

% LOOP OVER VARIOUS REPRESENTATIONS (COUPLINGS al,be)
ne = 4;                 % how many rep expts?
rhsnrm=cell(ne,1); munrm=cell(ne,1); gnrm=cell(ne,1); errt=cell(ne,1);
for r=1:ne, al = 1.0*mod(r-1,2); be = 2.0*(r>2);  % =====================
  fprintf('r=%d: al=%g, be=%g..........\n',r,al,be)
  als(r)=al; bes(r)=be;   % Representation is u = D.mu + be.S.mu + al.S.mudot:
  Rtarg = Dtarg + al*Sdottarg + be*Starg;  % for history application of ret BIEs
  Rtest = Dtest + al*Sdottest + be*Stest;  % for the test pt eval
  
  mo = []; mo.verb = 1; mo.shift=corrshift;  % opts for timestepping
  [tj u rhsnrm{r} gnrm{r} munrm{r} muall] = tmarch(dt,Ttot,predcorr,gdata,Rtarg,Rtest,wpred,mo);  % do time steps
  uex = data_ptsrc(xs,T,Tt,tj,t.x);            % known BVP soln at test pt
  maxerr(r) = max(abs(u-uex));
  errt{r} = u-uex;             % save error func vs t
  showsurffunc(s,muall(:,end)); title(sprintf('r=%d',r)); drawnow % show last mu
end                                              % ================== (r loop)

%clear Dtarg Starg Rtarg Sdottarg Stest Dtest Rtest Sdottest; save toggleab.mat
stop

chars = 'abcd'; % do all figs in one go................................
for r=1:ne
  figure; set(gcf,'position',[100+510*als(r),400-600*bes(r)/2,500,500]); 
  for s=1:2, subplot(2,1,s);
    plot(tj,abs(uex),'b-', tj,gnrm{r},'g:', tj,munrm{r},'r--');
    axis tight; v=axis; v(2)=27; axis(v);
    if s==1        % lin plot
      title(sprintf('(%s) coupling parameters $a=%g, b=%g$',chars(r),als(r),bes(r)),'interpreter', 'latex');
      if r==1, h=legend('$u$','$\|g(\cdot,t)\|_\infty$','$\|\mu(\cdot,t)\|_\infty$','location','northwest'); set(h,'interpreter','latex'); end
    else           % log plot
      hold on; plot(tj,abs(errt{r}),'k.-');    % add error (useless on lin plot)
      set(gca,'yscale','log'); v=axis; v(3:4)=[1e-12,1.1]; axis(v);
      yticks([1e-10,1e-5,1]);
      xlabel('t','interpreter','latex');
      if r==1, h=legend('$u$','$\|g(\cdot,t)\|_\infty$','$\|\mu(\cdot,t)\|_\infty$','$u$ error','location','southeast'); set(h,'interpreter','latex'); end
    end
  end
  set(gcf,'paperposition',[0 0 4 3.5]);
  print('-depsc2',sprintf('toggleab%d.eps',r));
end
  
surfarea =        22.6225981826694;   % for this shape (see end for eval)
fprintf('N=%d: h = %g\n',N,sqrt(surfarea/N))
%showsurffunc(s,muall(:,end));   % show unstable mode (if unstable)
