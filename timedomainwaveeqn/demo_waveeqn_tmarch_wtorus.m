% First baby timesteps in t-domain wave eqn BIE: known exterior Dirichlet BVP.
% Four different reps, output eps of max norms vs j, and videos of g, mu.
% Barnett 1/4/17-1/16/17, w/ Hagstrom, Greengard. 6/7/17 pred-corr.
% wobbly torus geom, 10/10/18.

% fsparse crashes for big probs (on desktop not laptop)
%rmpath ~/matlab/stenglib/Fast

%clear
dt = 0.1;   % timestep
m = 4;      % control time interp order (order actually m+2)
predcorr = -1;   % <0 for impl;  0,1,2, for pred-corr with that many corrector steps
wobbly = 1;
so.a=1; b=0.5; o.p=6;  % torus shape (a,b);  p = G-L nodes per panel side.
if ~wobbly  % plain torus
  so.b =b;
else
  wc = 0.1;  % surf modulation ampl
  wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
  wn = 5;   % # wobbles in toroidal, major
  f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
  ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
  so.b = {f,ft,fp};     % pass in instead of b param
end

% find dt=0.1 unstable for (6,4) pans of p=12, even tho (12,8) @ p=6 stable.
so.np=12; so.mp=so.np*2/3;    % panel numbers in major and minor direction.
%so.np=15; so.mp=10;
%so.np=6; so.mp=4;
% need to do a joint study over h and dt to check stability
[s N] = create_panels('torus',so,o); % surf
[x nx w] = getallnodes(s);
distmax = 6.0; %4.0;       % largest dist from anything to anything
n = ceil(distmax/dt);

if 1
SDload = 0;   % if load, params must match the above!
if SDload, disp('loading SD retarded history BIE matrices...')
  load SDtarg_wtorus_p6_m6_dt01      % precomputed 2GB (took 80 sec)
else
  o.nr = 2*o.p+4; o.nt = 2*o.nr;     % add aux quad to panels: aux quad orders
  % (note nr too small, eg 8, causes weak instability!)
  s = add_panels_auxquad(s,o);
  Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
  %profile clear; profile on
  t0=tic;
  [Starg,Dtarg,Sdottarg] = tdSDinterpmats_panels(s,s,Linfo,struct('n',n,'dt',dt,'m',m));
  fprintf('tot S,D,Sdot mat build (each %dx%d, nnz=%d): %.3g s\n',size(Starg,1),size(Starg,2),nnz(Starg),toc(t0))
  %profile off; profile viewer
  %disp('saving...'), save SDtarg_wtorus_p6_m6_dt01 Starg Dtarg Sdottarg -v7.3
end
end

% surf data for GRF test of SDtest vectors only...
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t); % data src t-func, tested
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside
% eval sig, tau on {n history grid} x {N bdry nodes}; recall time fast, N slow
tt = dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,T,Tt,ttt,xx,nxx);       % output ft unused
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

utest = Stest*sighist + Dtest*tauhist;   % two terms in GRF
uex = data_ptsrc(xs,T,Tt,0.0,t.x);     % exact u at test pt
fprintf('test that the u eval vectors (SDtest) work: %.3g\n', utest-uex)

mscs = [25 2 1.6 1.3];  % choice of max munows for color 3d plot

% LOOP OVER VARIOUS REPRESENTATIONS (SETTING al,be)
%for rep=1:4, al = mod(rep-1,2); be = rep>2; msc = mscs(rep);  %======= 1:4
for rep=4, al = 1.0*mod(rep-1,2); be = 2.0*(rep>2); msc = mscs(rep);  %======= 1:4, new al,be

% Representation is u = D.mu + be.S.mu + al.S.mudot:
%al = 1; be = 0;   % rep params, overridden by loop
Rtarg = Dtarg + al*Sdottarg + be*Starg;  % for history application of ret BIEs
Rtest = Dtest + al*Sdottest + be*Stest;  % for the test pt eval
%clear Dtarg Starg Sdottarg Dtest Stest Sdottest   % ? from now, R.. sparse mats

% pull out current-time matrix (last time index not first!):
Rnow = Rtarg(:,n:n:end); % NxN

if predcorr>=0   % set up for t-stepping
  wpred = extrap(m);          % row vector
  Cnow = diag(Rnow) + 1/2;    % col vector, note includes the 1/2
  Bnow = Rnow - diag(diag(Rnow));   % off-diag part
else        % implicit
  [Lnow Unow pnow] = lu(speye(N)/2 + Rnow,'vector'); % direct factor on-surf sys
%rhs = randn(N,1); mu = Unow\(Lnow\rhs(pnow)); norm(Rnow*mu + mu/2 - rhs) % check direct solve via LU
% [set up sparse to do extrap on each node separately (only for pred-corr)]
end
  
% Dirichlet data for BVP: t-dep pulse func for interior pt src... (max val ~1)
t0=6; s0=1.0; T = @(t) 5*exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
% (t0/s0 = 6 gives 1e-8 of start-up error if no time delay from src to surf)
%eps = 1e-5; tt = 4.3; fprintf('check Tt vs T: %.3g\n',(T(tt+eps)-T(tt-eps))/(2*eps) - Tt(tt)), clear tt

Ttot = 30.0;     % total time to evolve
jtot = ceil(Ttot/dt); 
verb = 1;    % 0: text only, 1: final plot, 2: anim, 3: save movie
muhist = zeros(n*N,1);  % init dens hist
gs=nan(jtot,1); rs=gs; ms=gs; es=gs;   % to save sizes of things for later
nam=sprintf('wtorus_p%d_m%d_dt01_pulse_al%gbe%g_march',o.p,m,al,be);
if verb>2, wO = VideoWriter([nam '.avi']); open(wO); end % mu scale
tic
for j=1:jtot           % .... marching loop
  tj = j*dt;
  muhist = reshape(muhist,[n N]); muhist = muhist([2:n,n],:); % shuffle back 1
  muhist(end,:)=0; muhist = reshape(muhist,[n*N 1]); % ensure munow=0
  gnow = data_ptsrc(xs,T,Tt,tj,x,nx);                % Dirichlet data from src
  rhs = gnow - Rtarg*muhist;                         % RHS for "now" lin sys
  if predcorr>=0
    muh = reshape(muhist,[n,N]); muh=muh(n-numel(wpred):n-1,:);
    munow = (wpred * muh)';  % predictor, kicks off iter for munow. row vec len N
    shift = 0.3;    % unknown best size
    t4=tic;
    for k=1:predcorr         % corrector steps, expressed via change in munow...
      dmunow = (rhs - Rnow*munow - munow/2) ./ (Cnow+shift);
      munow = munow + dmunow;
      fprintf('\t k=%d\t ||dmunow||=%g\t ||resid||=%g\n',k,norm(dmunow),norm(Rnow*munow+0.5*munow-rhs)) % ...so can track norm
      %munow = (rhs - Bnow*munow) ./ Cnow;  % plain corrector
    end
    fprintf('\t%d corr steps in %.3g s\n',predcorr,toc(t4))
  else    % implicit
    munow = Unow\(Lnow\rhs(pnow));           % write into "now" entries
  end
  %if mod(j,5)==0, showsurffunc(s,munow); title(fprintf('j=%d',j)), drawnow; end
  muhist(n:n:end) = munow;
  u = Rtest*muhist;                                  % eval at test pt
  uex = data_ptsrc(xs,T,Tt,tj,t.x);                  % known BVP soln at test pt
  fprintf('j=%d (tj=%g):\t err=%.3g \tu=%.6g \tuex=%.6g\n',j,tj,u-uex,u,uex)
  es(j)=u-uex; gs(j)=max(abs(gnow)); rs(j)=max(abs(rhs)); ms(j)=max(abs(muhist(n:n:end)));
  %imagesc(reshape(muhist,[n N])); caxis([0 3]); colorbar; drawnow
  if verb>1, %figure(1); plot([gnow, rhs, muhist(n:n:end)], '.-'); title(sprintf('t=%.3g',tj)); legend('g','rhs','\mu_{now}'); drawnow
    so=[]; so.nofig=1; figure(1); set(gcf,'position',[1000 500 500 800]); subplot(2,1,1);
    showsurffunc(s,gnow,so); light; caxis([0 1]); title(sprintf('g:  t=%4.1f',tj));
    subplot(2,1,2); showsurffunc(s,muhist(n:n:end),so); caxis([0 msc]); light; title('\mu_{now}'); drawnow;
    if verb>2, writeVideo(wO,getframe(gcf)); end
  end
end                    % ....
fprintf('done %d steps in %.3g sec: %.3g sec per t-step\n', jtot, toc, toc/jtot)
if verb>2, close(wO);     % writes AVI movie out; now encode small MP4...
  system(sprintf('ffmpeg -i %s.avi -y -c:v libx264 -crf 20 %s.mp4',nam,nam));
end
if verb, figure; semilogy([gs rs ms abs(es)],'.-');
  legend('||g||_\infty','||rhs||_\infty','||\mu_{now}||_\infty','u err','location','northwest');
  xlabel('timesteps'); axis tight;
  title(sprintf('wtorus, \\alpha=%g, \\beta=%g, non-osc pulse, \\delta t = %.3g, predcorr=%d',al,be,dt,predcorr));
  v=axis; v(3)=max(v(3),1e-11); axis(v); print('-depsc2',[nam '_log.eps'])
  %set(gca,'yscale','lin'); print('-depsc2',[nam '.eps'])
end

end % ================== (rep loop)


% show unstable mode (if unstable)...
%showsurffunc(s,muhist(n:n:end));
