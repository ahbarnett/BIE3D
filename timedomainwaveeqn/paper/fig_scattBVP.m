% t-domain wave eqn BIE: scattering (exterior Dirichlet BVP).
% paper fig version, Barnett 1/9/19, based on fig_toggleab.m

% fsparse crashes for big probs (on desktop not laptop)
rmpath ~/matlab/stenglib/Fast
clear; verb = 3;

dt =   0.05;   % timestep
Ttot = 10.0;     % total time to evolve
m = 4;      % control time interp order (order actually m+2)
predcorr = 8;   % <0 for impl;  0,1,2, for pred-corr with that many corr steps
corrshift = 0.25;    % only helps for large dt
wobbly = 1;   % 0: torus, 1: cruller.
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
so.np = 18; so.mp=round(so.np*2/3);  % panel # in major,minor directions
o.nr = 2*o.p; o.nt = 2*o.nr;  % aux quad orders: radial and angular nodes
% Dirichlet data for BVP: t-dep pulse func...
incwave = 'plane';    % 'ptsrc' (lauched at t0) or 'plane' (hitting 0 at t0)
s0 = 0.12;         % closely tied to dt. Error ~ exp(-3*(dt/s0)^2)
if strcmp(incwave,'ptsrc')
  a0=50; t0=-1; [T,Tt] = pulsegaussian(a0,t0,s0);
  xs = [.7;0.3;3.5];   % src pt for data, above
  uinc = @(x,t) data_ptsrc(xs,T,Tt,t,x);    % incident wave func
elseif strcmp(incwave,'plane')
  a0=1; t0=3; [T,Tt] = pulsegaussian(a0,t0,s0);
  incdir = [-.2,.1,-1]; incdir = incdir/norm(incdir);
  uinc = @(x,t) data_planewave(incdir,T,Tt,t,x);    % incident wave func
end
fwhm = fzero(@(t) T(t)-a0/2,t0+1)-fzero(@(t) T(t)-a0/2,t0-1)  % get pulse width

shape='torus'; [s N] = create_panels(shape,so,o);     % surf (torus class)
surfarea =        22.6225981826694;   % for this shape
h = sqrt(surfarea/N);
[x nx w] = getallnodes(s);
gdata = @(t) -uinc(x,t);    % cancel uinc Dirichlet data on s.x, as func of t
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
ttest=6.0; tt = ttest + dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,T,Tt,ttt,xx,nxx);       % (output ft unused)
sighist = -fn; tauhist = f;  % col vecs, ext wave eqn GRF: u = D.u - S.un
%f = reshape(f,[n N]); figure;
%oo.nofig=1; oo.sc=1; for i=1:n, showsurffunc(s,f(i,:),oo); drawnow, end

% set up vectors which compute potential u(t,x) at >1 targs x, from dens hist...
%t.N = 100; t.x = [1.0;0.2;0.5] + [0;0;3]*(1:t.N)/t.N;  % line of targs
dx=0.05; gx=-2:dx:2; gz = -1:dx:2; [xx zz] = meshgrid(gx,gz);
t.x = [xx(:)';0*xx(:)';zz(:)']; t.N=numel(xx);           % x-z plane of targs
showsurffunc(s,0*x(1,:)); hold on; plot3(xs(1),xs(2),xs(3),'r.');
plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.'); drawnow

t0=tic;
tret = -dists(t.x,x);  % retarded times of surf nodes rel to ttarg, trg loc fast
[jmax,jmin,a,ap] = interpmat(tret,dt,m);   % Tom's coeffs, 1 row per tret entry
clear tret;   % hitting 30 GB even at np=9, why?
joff = jmin+n-1;         % padding on the ancient side
if joff<0, error('u_test eval vec requesting too ancient history!'); end
a = [zeros(N*t.N,joff), a, zeros(N*t.N,-jmax)];  % dense, pad to width n
ap = [zeros(N*t.N,joff), ap, zeros(N*t.N,-jmax)];
[S D Dp] = tdSDmats(t.x,x,nx,w);  % each is NtxN (Nt = N.t = # targ)
Stest = a.*repmat(S(:),[1 n]);  % coeff vectors packed as (Nt*N) x n mats
Stest = sparse(reshape(permute(reshape(Stest,[t.N N n]),[1 3 2]),[t.N N*n]));
Dtest = a.*repmat(D(:),[1 n]) + ap.*repmat(Dp(:),[1 n]); clear a;
Dtest = sparse(reshape(permute(reshape(Dtest,[t.N N n]),[1 3 2]),[t.N N*n]));
Sdottest = ap.*repmat(S(:),[1 n]); clear ap;
Sdottest = sparse(reshape(permute(reshape(Sdottest,[t.N N n]),[1 3 2]),[t.N N*n]));
fprintf('test target eval mats filled, %.3g s\n',toc(t0));
% this test only makes sense if targs far from surf...
utest = Stest*sighist + Dtest*tauhist;   % two terms in eval GRF
fprintf('GRF null-field test, max at all %d targs: %.3g\n', t.N, max(abs(utest)))

al=1; be=2;   % Representation is u = D.mu + be.S.mu + al.S.mudot:
Rtarg = Dtarg + al*Sdottarg + be*Starg;   % for history application of ret BIEs
Rtest = Dtest + al*Sdottest + be*Stest;   % (sparse) for the test pt eval
clear Stest Dtest Sdottest Starg Dtarg Sdottarg;
% EVOLVE --------
mo = []; mo.verb = 1; mo.shift=corrshift;  % opts for timestepping
[tj u rhsnrm gnrm munrm muall] = tmarch(dt,Ttot,predcorr,gdata,Rtarg,Rtest,wpred,mo);  % do time steps; u is now scattered BVP soln
nst = numel(tj);    % # time steps done
ttt = kron(tj',ones(1,t.N)); xxx = kron(ones(1,nst),t.x); % spacetime eval lists
utot = reshape(uinc(xxx,ttt) + u(:), [t.N nst]);      % #targs * #tsteps
% ---------------

% check at point...
t0=4.0; jt=find(tj==t0); jx=37; jz=41; j=jz+numel(gz)*(jx-1);
fprintf('test pt x0=(%g,%g,%g)\n',xx(j),0,zz(j));
fprintf('dt=%g,m=%d,np=%d,p=%d,nr=%d: \t u(x0,%g) = %.9f\n',dt,m,so.np,o.p,o.nr,t0,utot(j,jt))
mask = ~insidexsecdomain(shape,so,xx,zz);   % masks out interior pts in slice

if verb>1   % rect xz-slice (y=0) array of targets, anim
  figure; for i=1:nst
  imagesc(gx,gz,mask.*reshape(utot(:,i),size(zz))); xlabel('x'); ylabel('z');
  caxis(2*[-1 1]); axis xy equal tight; v=axis; h=showsurfxsec(shape,so);
  axis(v); title(sprintf('t=%.2f',tj(i))); drawnow; pause(0.05); hold off;
  end
  figure; plot(tj,utot(jz+numel(gz)*(jx-1),:),'+-'); xlabel('t');
end
  
if 0 % line target history...
figure; surf(tj,t.x(3,:),utot); shading interp; colorbar; axis vis3d
ylabel('z coord of targ'); xlabel('t time'); v=axis; v(3)=0.5;
figure; for i=1:nst              % anim for line targets
  plot(t.x(3,:),utot(:,i),'-'); sc=max(abs(utot(:)));
  axis([.5 max(t.x(3,:)) -sc sc]); title(sprintf('t=%.3f',tj(i)));
drawnow; pause(0.02); end
end

if verb>1  % 3d slice anim...
  nam=sprintf('cruller_scatt_%s_pulse_dt%g_m%d_p%d_np%d',incwave,dt,m,o.p,so.np);
  if verb>2, wO=VideoWriter([nam '.avi']); wO.FrameRate=15; wO.open; end
  figure; for i=1:nst
  oo=[]; oo.nofig=1; [h0 h4]=showsurffunc(s,muall(:,i),oo); hold on;
  h1=surf(xx,0*xx,zz,mask.*reshape(utot(:,i),size(zz)));
  set(h1,'FaceLighting','none','linestyle','none','facecolor','flat');
  oo=[]; oo.dims=3; h2=showsurfxsec(shape,so,oo);   % add intersection curves
  caxis(1.0*[-1 1]); view(-30+i/2,35); lightangle(45,0); axis tight;
  ax=gca; ax.Clipping='off'; ax.CameraPosition = ax.CameraPosition*0.6;
  h3=title(sprintf('slice of u_{tot} and \\mu on surface: t=%.2f',tj(i)));
  h3.Position = [0,0,2.5];
  pos = h4.Position; h4.Position = pos + [.07,0,0,0];  % push colorbar right
  hold off; drawnow;
  if verb>2, writeVideo(wO,getframe(gcf)); end
  if i<nst, clf; end   % otherwise surf not cleared, unsure why
  end
  if verb>2, close(wO);     % writes AVI movie out; now encode small MP4...
    system(sprintf('unset LD_LIBRARY_PATH; ffmpeg -i %s.avi -y -c:v libx264 -crf 20 %s.mp4',nam,nam));
  end
end

% ============================================================================
if verb   % figs for paper... (will have to output PNG then convert to EPS)
  jx=37; jz=41; j=jz+numel(gz)*(jx-1); x0=xx(j); z0=zz(j);  % spatial pt
  fprintf('test pt x0=(%g,%g,%g)\n',x0,0,z0);
  t0s = [3.0, 7.0];           % two times, for (a) and (b)
  for ii=1:2
    figure; t0=t0s(ii);
    i = find(abs(tj-t0)<1e-14);
    oo=[]; oo.nofig=1; [h0 h4]=showsurffunc(s,muall(:,i),oo); hold on;
    h1=surf(xx,0*xx,zz,mask.*reshape(utot(:,i),size(zz)));
    set(h1,'FaceLighting','none','linestyle','none','facecolor','flat');
    oo=[]; oo.dims=3; h2=showsurfxsec(shape,so,oo);   % add intersection curves
    plot3(x0,0,z0,'.','markersize',20);
    caxis(1.0*[-1 1]); if ii==2, caxis(max(abs(utot(:,i)))*[-1 1]); end
    view(20,35); lightangle(45,0); axis tight;
    ax=gca; ax.Clipping='off'; ax.CameraPosition = ax.CameraPosition*0.65;
    h3=title(sprintf('(%c)     u_{tot} on \\{y=0\\} and \\mu on \\Gamma:     t=%.g',ii+96,tj(i)));
    h3.Position = [0,0,2.5];
    pos = h4.Position; h4.Position = pos + [.07,0,0,0];  % push colorbar right
    set(gcf,'paperposition',[0 0 5.8 5]);
    nam = sprintf('scatt_%d',ii);
    print('-dpng','-r600', [nam '.png']);
    system(['convert -trim ' nam '.png eps2:' nam '.eps']);
  end
  l = load('../expts/wobblytorus/scattBVPconv.mat');  % by:gen_scattBVPconv.m
  r = l.run{3};
  figure; plot(r.tj,r.utot(jz+numel(gz)*(jx-1),:),'.-');   % signal
  xlabel('$$t$$','interpreter','latex');
  ylabel('$$u_{tot}(x_0,t)$$','interpreter','latex');
  title('(c) \quad total wave signal at $$x_0$$', 'interpreter','latex');
  set(gcf,'paperposition',[0 0 5 3]);
  print -depsc2 scatt_sig.eps
  
  % (d) error convergence at x_0, t=3,  vs np = 6,9
  % with dt = h.
  % **** loop over runs
end %=======
