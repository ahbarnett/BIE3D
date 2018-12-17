% Explore stability in dx-dt plane for explicit t-step wave eqn BIE.
% wobblytorus shape, interior source RHS (known solution) for now.
% code inside loop from waveeqn_tmarch.m, now calls tmarch.m
% Dependencies: rest of BIE3D, memorygraph, optional:Parallel Toolbox

run('../../../bie3dsetup.m')   % run from expts/wobblytorus
% either use or skip fsparse (crashed on ccblin019 for np>=15)...
%addpath ~/matlab/stenglib/Fast      % danger: only use for small problems!
rmpath ~/matlab/stenglib/Fast
addpath ~/matlab/memorygraph

clear all   % needed since tstepmatvec has persistent, shoul
so.a=1; b=0.5; % torus shape (a,b)
wc = 0.1;  % surf modulation ampl (0 for plain torus)
wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
wn = 5;   % # wobbles in toroidal, major
f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
so.b = {f,ft,fp};     % pass in instead of b param
distmax = 6.0; %4.0;       % largest dist from anything to anything (diam).

% Dirichlet data for BVP: t-dep pulse func for interior pt src... (max val ~1)
xs = [0.9;-0.2;0.1];   % src pt for RHS, must be inside
t0=6; s0=1.0; T = @(t) 5*exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
% (t0/s0 = 6 gives 1e-8 of start-up error if no time delay from src to surf)
%eps = 1e-5; tt = 4.3; fprintf('check Tt vs T: %.3g\n',(T(tt+eps)-T(tt-eps))/(2*eps) - Tt(tt)), clear tt

Ttot = 50.0; %16; % total time to evolve (pulse has passed by 16; want tail)

o.p = 8;  % p = spatial "order" (G-L nodes per panel side)
m = 6;    % control time interp order (order actually m+2)
wpred = extrap(m);   % (m+2)th order extrapolation row vector for prediction

al = 1.0;    % repr mixing alpha
nam = sprintf('wt_T50_m%d_p%d_pulse',m,o.p);

% define dx and dt parameter plane...
nps = 6:3:18;   % num panels on major torus loop (multiples of 3 best)
%dts = [0.03 0.04 0.05 0.06 0.07 0.085 0.1 0.12 0.15 0.2 0.3 0.4 0.5];  % dt timesteps
%dts = 0.1*2.^(-2:0.5:3);
dts = 0.1*2.^(-2:0.25:3);
maxsteps=max(ceil(Ttot./dts));

% define expt runs for each dx-dt pair...
betas = 2; %[0.5 1 2 4];    % beta repr mixing params
pcs = [8]; %[1,2,3,4,-1];     % pred-corr steps (1,2,..;-1=inf)
corrshift = 0.25;
nbe = numel(betas); npc = numel(pcs);

ne = nbe*npc;  % # expts for each (dx,dt) pair, ie S,D full matrices filled
gros = nan(numel(dts),numel(nps),nbe,npc); % where measured growth factors will go
errs = gros;      % sup errors
errts = nan(maxsteps,numel(dts),numel(nps),nbe,npc);  % 5d t-steps output arrs
gnrms=errts; rhsnrms=errts; munrms=errts;
neigful = 3;                      % # eigs of full t-step matrix
neigjac = 12;                  % to save spectrum of jacobi matrix
eigful = nan(neigful,numel(dts),numel(nps),nbe,npc);  % full t-step mat
eigmax = nan(numel(dts),numel(nps),nbe,npc);  % est max abs eig
eigjac = nan(numel(dts),numel(nps),2,neigjac);
memorygraph('start');

for i=1:numel(nps)    % ======================= MAIN DX LOOP

  t3=tic;
  so.np=nps(i); so.mp = round(so.np/3*2);   % spatial discr panel numbers
  [s N] = create_panels('torus',so,o); % surf
  [x nx w] = getallnodes(s);
  o.nr = 12; o.nt = 2*o.nr;    % first add aux quad to panels: aux quad orders
  s = add_panels_auxquad(s,o);
  Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad

  for j=1:numel(dts), dt=dts(j);    % --------------------- MAIN DT LOOP
    n = ceil(distmax/dt);   % # history for apply ops
    fprintf('p=%d (np,mp)=(%d,%d) N=%d\tdt=%.3g n=%d --------\n',o.p,so.np,so.mp,N,dt,n)
    tt=tic;  % expensive setup space-time history self-quadrature on surface...
    [Starg,Dtarg,Sdottarg] = tdSDinterpmats_panels(s,s,Linfo,struct('n',n,'dt',dt,'m',m));  % NB uses parfor
    fprintf('\tS,D,Sdot ret mat build (each %dx%d, nnz=%d): %.3g s\n',size(Starg,1),size(Starg,2),nnz(Starg),toc(tt))

    % set up vectors which compute potential at u(tnow,x_ext) from dens hist...
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
    % (see waveeqn_tmarch for GRF test of these test-pt-eval vecs)

    errsji = nan(1,ne); grosji = nan(1,ne); % must be 1d for parfor
    rhsnrm=cell(ne,1); munrm=cell(ne,1); gnrm=cell(ne,1); errt=cell(ne,1);
    eigfulji=cell(ne,1); eigmaxji = nan(1,ne);
    t2=tic; fprintf('----- get ready for %d timestepping expts -----\n',ne)
    for e=1:ne         % if parfor, restricts outputs to be 1d in e
      ib = 1+mod(e-1,nbe); ip = 1+floor((e-1)/nbe);  % decode e into 2 inds
      be = betas(ib);
      Rtarg = Dtarg + al*Sdottarg + be*Starg;  % for history appl of ret BIEs
      Rtest = Dtest + al*Sdottest + be*Stest;  % for the test pt eval
      if e==1      % for Rnow solve, output 2-by-neigjac (assumes beta fixed)
        eigjacji = specjacobi(Rtarg(:,n:n:end),corrshift,neigjac);
      end
      predcorr = pcs(ip);
      fprintf('expt e=%d:\t beta=%.3g \t predcorr=%d ..........\n',e,be,predcorr)
      gdata = @(t) data_ptsrc(xs,T,Tt,t,x,nx);  % Dirichlet data from xs src
      mo = []; mo.verb = 1; mo.shift=corrshift;  % opts for timestepping
      
      if 0      % max abs eigval estimation (L = exp(gro)),  2 ways...
        matvec = @(muhist) tstepmatvec(muhist, predcorr, Rtarg, wpred, mo);
        clear tstepmatvec      % resets this func's persistent vars
        t5=tic;
        disp('arpack eigs of tstep mat...');   % ~20 matvecs per iter...
        [V,D,flag] = eigs(matvec,n*N,neigful,'LM',struct('tol',1e-4,'maxit',50,'disp',0));
        fprintf('arpack eigs done in %.3g s (flag=%d): largest eig mag = %.6g\n',toc(t5),flag,max(abs(diag(D))));
        eigfulji{e} = diag(D);
        if max(abs(diag(D)))<1.1    % ambiguous
          clear tstepmatvec      % resets this func's persistent vars
          t5=tic;
          disp('power method on tstep mat...');
          eigmaxji(e) = maxeigpowermeth(matvec,n*N,500,100);  % enough
          fprintf('maxeigpowermeth done in %.3g s: largest eig mag = %.6g\n',toc(t5),eigmaxji(e));
        end
      end
      
      %mo.random = 1;   % estimate largest eigenval by power method (old)
      [tj u rhsnrm{e} gnrm{e} munrm{e}] = tmarch(dt,Ttot,predcorr,gdata,Rtarg,Rtest,wpred,mo);  % time step
      uex = data_ptsrc(xs,T,Tt,tj,t.x);            % known BVP soln at test pt
      errsji(e) = max(abs(u-uex));
      errt{e} = u-uex;             % save error func vs t
      pad=20;                     % how many dt steps to use to meas mu gro rate
      grosji(e) = log(munrm{e}(end)/munrm{e}(end-pad))/pad; % rate~e^(grow.t/dt)
      fprintf('     e=%d found err=%.3g \t gro=%.3g\n',e,errsji(e),grosji(e))
    end
    fprintf('----- %d timestepping expts done in %.3g s total -----\n',ne,toc(t2))
    % extract 1d stuff from parfor into full output arrays...
    errs(j,i,:,:) = reshape(errsji,nbe,npc);
    gros(j,i,:,:) = reshape(grosji,nbe,npc);
    eigjac(j,i,:,:) = eigjacji;
    for e=1:ne, ib = 1+mod(e-1,nbe); ip = 1+floor((e-1)/nbe);   % decode e
      nt = ceil(Ttot/dt); rhsnrms(1:nt,j,i,ib,ip) = rhsnrm{e};
      errts(1:nt,j,i,ib,ip) = errt{e};
      munrms(1:nt,j,i,ib,ip) = munrm{e};
      gnrms(1:nt,j,i,ib,ip) = gnrm{e};
      %eigful(1:neigful,j,i,ib,ip) = eigfulji{e};
      eigmax(j,i,ib,ip) = eigmaxji(e);
    end
    clear errt munrm gnrm rhsnrm Starg Dtarg Sdottarg Rtarg   % kill big mats before rebuild them
  
    if 0 % example plots: all errors (for this dt, and iff single beta).
         %j = 1;  % choose the dt to study
      figure; semilogy(tj,abs(squeeze(errts(1:nt,j,i,1,:))),'+-');
      legnum(pcs); hold on; set(gca,'ColorOrderIndex',1)
      semilogy(tj,squeeze(munrms(1:nt,j,i,1,:)),'.--');
      axis tight; v=axis; v(3:4)=[1e-20,1];
      semilogy(tj,gnrms(1:nt,j,i,1,1),'k:'); axis(v);
      title(sprintf('wtorus p=%d np=%d m=%d dt=%g',o.p,nps,m,dt))
    end
  end                               % ---------------------

  tnps(i) = toc(t3);
  fprintf('entire np=%d case done in %.3g s\n',nps(i),tnps(i))
  % whos  % useful for log files
  % now kill big stuff before save...
  clear Starg Dtarg Sdottarg Rtarg Dtest Stest Sdottest Rtest a ap s Linfo
  memorygraph('label',sprintf('np=%d done',nps(i)));
  [bytes est_times cpu_times cpu_usages labelstrings labeltimes] = memorygraph('get');
  save(nam)   % save after every dx surf discr choice
  
end                   % =======================

memorygraph('done');
%memorygraph('plot',bytes,est_times,cpu_times,cpu_usages,labelstrings,labeltimes);

%return   % needed in batch mode? no, quit. But if use matlab < script.m,  ok
% don't use matlab -r script.


% test of stability  -------------------------------------------------

if 0 % compare the 3 ways of estimating max eigval of H, ie (in)stability...
  1+squeeze(gros)'
  max(squeeze(abs(eigful(:,1,1,1,:))))
  squeeze(eigmax(1,1,1,:))'
  max(squeeze(abs(eigful(:,1,1,1,:)))) - squeeze(eigmax(1,1,1,:))'
end
  
if 0 % spectra (color = #pc) from arpack...
  figure; plot(squeeze(eigful(:,1,1,1,:)),'+'); axis equal; hold on;
  plot(exp(2i*pi*(1:1e3)/1e3),'m-');
end
