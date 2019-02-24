% t-domain wave eqn BIE: scattering (exterior Dirichlet BVP). Convergence.
% paper fig version, Barnett 1/14/19.
% dx coupled with dt, convergence at point (x0,0,z0).

run('../../../bie3dsetup');
% fsparse crashes for big probs (on desktop not laptop)
rmpath ~/matlab/stenglib/Fast
clear;

Ttot = 10.0;     % total time to evolve
m = 4;      % control time interp order (order actually m+2)
predcorr = 8;   % <0 for impl;  0,1,2, for pred-corr with that many corr steps
corrshift = 0.25;    % only helps for large dt
shape='torus';         % shape class
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
% Dirichlet data for BVP: t-dep pulse func...
incwave = 'plane';    % 'ptsrc' (lauched at t0) or 'plane' (hitting 0 at t0)
s0 = 0.12;         % closely tied to dt. Error ~ exp(-3*(dt/s0)^2)
if strcmp(incwave,'ptsrc')
  a0=50; t0=-1; [T,Tt] = pulsegaussian(a0,t0,s0);
  xs = [.7;0.3;3.5];   % src pt for data, above
  uinc = @(x,t) data_ptsrc(xs,T,Tt,t,x);    % incident wave func
elseif strcmp(incwave,'plane')
  a0=1; t0=3; [T,Tt] = pulsegaussian(a0,t0,s0);
  incdir = [-.2,.1,-1]; incdir = incdir/norm(incdir);  % shine wave down
  uinc = @(x,t) data_planewave(incdir,T,Tt,t,x);    % incident wave func
end
fwhm = fzero(@(t) T(t)-a0/2,t0+1)-fzero(@(t) T(t)-a0/2,t0-1)  % get pulse width

distmax = 6.0;                   % largest dist from anything to anything
surfarea =        22.6225981826694;   % for this shape
o.nr = 2*o.p; o.nt = 2*o.nr;  % aux quad orders: radial and angular nodes

nps = 6:3:24;    % convergence

memorygraph('start');
for i=1:numel(nps)         % ======================== MAIN DX LOOP ========
  t3 = tic;
  so.np = nps(i); so.mp=round(so.np*2/3);  % panel # in major,minor directions

  [s N] = create_panels(shape,so,o);     % surf (torus class)
  h = sqrt(surfarea/N);             % typical dx
  [x nx w] = getallnodes(s);
  gdata = @(t) -uinc(x,t);    % cancel uinc Dirichlet data on s.x, as func of t
  wpred = extrap(m);          % extrapolation row vector

  dt = 0.9 * h;      % choose timestep!  (note want small to capture pulse)
  % *** be ready to adjust down a bit, to 0.8 (close to inv CFL)

  fprintf('run %d: p=%d, %dx%d panels, N=%d, dt=%.3g\n',i,o.p,so.np,so.mp,N,dt)
  n = ceil(distmax/dt);
  s = add_panels_auxquad(s,o);
  Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
  t0=tic;
  [Starg,Dtarg,Sdottarg] = tdSDinterpmats_panels(s,s,Linfo,struct('n',n,'dt',dt,'m',m));
  fprintf('tot S,D,Sdot mat build (each %dx%d, nnz=%d): %.3g s\n',size(Starg,1),size(Starg,2),nnz(Starg),toc(t0))

  % set up vecs which compute potential u(t,x) at >1 targs x, from dens hist..
  dx=0.05; gx=-2:dx:2; gz = -1:dx:2; [xx zz] = meshgrid(gx,gz);
  t.x = [xx(:)';0*xx(:)';zz(:)']; t.N=numel(xx);           % x-z plane of targs
  %showsurffunc(s,0*x(1,:)); hold on; plot3(xs(1),xs(2),xs(3),'r.');
  %plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.'); drawnow

  t0=tic;
  if 0 % old, all mu-hist to test targets in a single dense op, scales poorly...
  tret = -dists(t.x,x); % ret times of surf nodes rel to ttarg, trg loc fast
  [jmax,jmin,a,ap] = interpmat(tret,dt,m); % Tom's coeffs, 1 row per tret entry
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
  else  % lower-RAM sparse fill of Stest, row by row, single assemble, etc...
    nnzm = ceil(1.1*(m+2)*t.N*N);  ii1 = zeros(1,nnzm);    % sparse index lists
    jj1=ii1;v1=ii1;ii2=ii1;jj2=ii1;v2=ii1;ii3=ii1;jj3=ii1;v3=ii1;
    p1=0;p2=0;p3=0;   % ptrs to sparse ind lists
    for it=1:t.N        % loop over test targs
      tret = -dists(t.x(:,it),x);   % ret times of surf nodes rel to ttarg
      [jmax,jmin,a,ap] = interpmat(tret,dt,m);  % Tom's, 1 row per tret entry
      joff = jmin+n-1;         % padding on the ancient side
      if joff<0, error('u_test eval vec requesting too ancient history!'); end
      a = sparse(a); ap = sparse(ap);
      a = [zeros(N,joff), a, zeros(N,-jmax)];   % sparse, pad to N*n
      ap = [zeros(N,joff), ap, zeros(N,-jmax)];
      [S D Dp] = tdSDmats(t.x(:,it),x,nx,w);    % each is 1xN, dense row
      z = a'.*kron(ones(n,1),S);  % a is sparse, S dense expanded vertically
      [~,jj,vv] = find(z(:)'); ind = p1+(1:numel(jj));  % unroll n*N sparse row
      ii1(ind)=it; jj1(ind)=jj; v1(ind)=vv; p1=ind(end);  % append
      z = a'.*kron(ones(n,1),D) + ap'.*kron(ones(n,1),Dp);
      [~,jj,vv] = find(z(:)'); ind = p2+(1:numel(jj));  % unroll n*N sparse row
      ii2(ind)=it; jj2(ind)=jj; v2(ind)=vv; p2=ind(end);  % append
      z = ap'.*kron(ones(n,1),S);
      [~,jj,vv] = find(z(:)'); ind = p3+(1:numel(jj));  % unroll n*N sparse row
      ii3(ind)=it; jj3(ind)=jj; v3(ind)=vv; p3=ind(end);  % append
    end
    Stest = sparse(ii1(1:p1),jj1(1:p1),v1(1:p1),t.N,N*n); % assemble in one step
    clear ii1 jj1 v1 p1 tret jj vv a ap ind
    Dtest = sparse(ii2(1:p2),jj2(1:p2),v2(1:p2),t.N,N*n);
    clear ii2 jj2 v2 p2
    Sdottest = sparse(ii3(1:p3),jj3(1:p3),v3(1:p3),t.N,N*n);
    clear ii3 jj3 v3 p3
  end
  fprintf('test target eval mats filled, %.3g s\n',toc(t0));

  al=1; be=2;   % Representation is u = D.mu + be.S.mu + al.S.mudot:
  Rtarg = Dtarg + al*Sdottarg + be*Starg;   % for hist application of ret BIEs
  Rtest = Dtest + al*Sdottest + be*Stest;   % (sparse) for the test pt eval
  clear Stest Dtest Sdottest Starg Dtarg Sdottarg;
  % EVOLVE --------
  mo = []; mo.verb = 1; mo.shift=corrshift;  % opts for timestepping
  [tj u rhsnrm gnrm munrm muall] = tmarch(dt,Ttot,predcorr,gdata,Rtarg,Rtest,wpred,mo);  % do time steps; u is now scattered BVP soln
  nst = numel(tj);    % # time steps done
  ttt = kron(tj',ones(1,t.N)); xxx = kron(ones(1,nst),t.x); % spacetime eval lists
  utot = reshape(uinc(xxx,ttt) + u(:), [t.N nst]);      % #targs * #tsteps
  % ---------------

  % check at point... (x0,0,z0) and t0
  t0=3.5; jt=find(abs(tj-t0)==min(abs(tj-t0)));  % rough similar pt
  jx=37; jz=41; j=jz+numel(gz)*(jx-1); x0=xx(j); z0=zz(j);  % spatial pt
  fprintf('test pt x0=(%g,%g,%g)\n',x0,0,z0);
  fprintf('dt=%g,m=%d,np=%d,p=%d,nr=%d: \t u(x0,%g) = %.9f (not conv test!)\n',dt,m,so.np,o.p,o.nr,tj(jt),utot(j,jt))   % note since dt) is off-dt-grid.
  mask = ~insidexsecdomain(shape,so,xx,zz);   % masks out interior pts in slice

  tnps(i) = toc(t3);
  fprintf('entire np=%d case done in %.3g s\n',nps(i),tnps(i))
  % whos  % useful for log files
  % now kill big stuff before save...
  clear Rtarg Rtest s Linfo D Dp S xxx ttt
  memorygraph('label',sprintf('np=%d done',nps(i)));
  [bytes est_times cpu_times cpu_usages labelstrings labeltimes] = memorygraph('get');

  run{i}.utot = utot; run{i}.tj = tj; run{i}.t = t;  % keep some stuff
  run{i}.j = j;  % space index of utot to look at, via t-interp to t0  ***
  run{i}.t0 =t0; run{i}.muall = muall; run{i}.dt = dt; run{i}.h = h;
  run{i}.np = nps(i);
  run{i}.maxbytes = max(bytes);

  save scattBVPconv.mat   % save after every dx surf discr choice
  
end                   % =======================

memorygraph('done');
