% spatial quadr scheme convergence, for paper. Barnett 11/22/18-12/12/18
% is most of test_tdGRF.m with convergence loop: it tests the retarded GRF.

clear; wobbly = 1;
so.a = 1; b = 0.5;   % major radius and poloidal radius
if ~wobbly  % plain torus
  so.b =b;
else         % cruller
  wc = 0.1;  % surf modulation ampl
  wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
  wn = 5;   % # wobbles in toroidal, major
  f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
  ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
  so.b = {f,ft,fp};     % pass in instead of b param
end
surfarea =        22.6225981826694;   % for this shape (see end for eval)

w0 = 5.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t);  % time data driving
Tfunctest = 0; if Tfunctest
  eps = 1e-5; t = 0.3;       % test time
  fprintf('error in time deriv: %.3g\n',(T(t+eps)-T(t-eps))/(2*eps) - Tt(t))
  clear t
end
xs = [0.9;-0.2;0.1];   % src pt for data, must be well inside
xext = [1.3;0.1;0.8];   % exterior targ
t.N = 1; ttarg = 2.1;   % test target time (avoid s.t conflict)

for fig=1:2, fig     %%%%%%% subfigure loop over p
  % (old figs had nps up to 27 for p=4, 18 for p=8...)
  if fig==1, o.p=4; nps = 6:3:33; nrs = [nan 8 12 16]; % order, np conv, nr conv
  else,      o.p=8; nps = 6:3:24; nrs = [nan 16 20 24];   % "
  end
  sides = [1 0 0 0];  % to match nr cases, must always be four
  errs = nan(numel(nps),numel(nrs));
  for c=1:numel(nrs)   % ======= lines (cases) on same plot, nr loop =========
    nr = nrs(c), side = sides(c);
    %side = 0;  % -1,0,1: choose nature of GRF test pt (will be a fake 1-pt panel)
    % (side=1 cleaner convergence, but for nr=12 or 14, side=0 is pretty close)

    Ns = nan*nps; SAs = Ns;
    for e = 1:numel(nps)         % ---------- np convergence -----------------
      so.np = nps(e); so.mp = round(2/3*so.np); so.orig = 0;  % shifts to ph=th=0
      [s N] = create_panels('torus',so,o); Ns(e)=N;
      if side==1, t.x = xext; nam='exterior';
      elseif side==-1, t.x = [-.8;-.3;.2]; nam='interior';
      else, %k=numel(s);j=o.p^2/2-o.p/2;   % case so.orig=1 (default offset)
        k=1; j=1; %o.p^2/2-o.p/2;   % node close to phi=th=0 (if so.orig=0), extreme tri's
        t.x = s{k}.x(:,j); % ...or on-surf (trg moves w/ np!)
      end
      %t.x  % check the targ is roughly const as converge

      [x nx w] = getallnodes(s);
      tret = ttarg - dists(t.x,x);   % retarded times on surf rel to GRF test pt
      if side~=0          % ------------------- off-surf 
        % now eval retarded sig, tau, tau'...
        [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);
        retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
        % do u(x,t) = S.sigma + (wave eqn D).tau :
        u = LapSeval_panels(t,s,retsig) + LapDeval_panels(t,s,rettau) + RetDeval_panels(t,s,rettaut);
        uex = data_ptsrc(xs,T,Tt,ttarg,t.x) * (side==1);      % zero inside, f out
        fprintf('N=%d, %s GRF test at 1 pt: u err = %.3g\n', N, nam, u-uex)
        SAs(e) = sum(w);
        
      else  % ------------- on-surf GRF w/ local aux nodes, target pan k, pt j
        o.nr = nr; o.nt = 2*o.nr;
        s = add_panels_auxquad(s,o);    % setup aux quad
        s = add_panels_timequad(s);     % local aux t-delays
        kill = [k; s{k}.nei];     % pan inds to remove from far list
        sf = s; sf(kill) = [];    % far pans - need fix internal inds in sf ???
        [x nx w] = getallnodes(sf);
        tret = ttarg - dists(t.x,x);   % retarded times on surf rel to GRF test pt
        [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);      % far
        retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
        u = LapSeval_panels(t,sf,retsig) + LapDeval_panels(t,sf,rettau) + RetDeval_panels(t,sf,rettaut);
        % make fake near panel w/ all aux nodes and hack the time delays for now..
        sn.x = s{k}.auxnodes(:,:,j);  sn.nx = s{k}.auxnormals(:,:,j);
        sn.w = s{k}.auxwei(:,j)';           % must be col vec
        sn.N = numel(sn.w);       % finish the fake aux src pan
        tret = ttarg - s{k}.auxdelays(:,j)';   % near panel's aux eval t's col vec
        [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,sn.x,sn.nx);      % near
        retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
        u = u + LapSeval_panels(t,sn,retsig) + LapDeval_panels(t,sn,rettau) + RetDeval_panels(t,sn,rettaut);
        uex = data_ptsrc(xs,T,Tt,ttarg,t.x) / 2;       % on-surf GRF
        fprintf('N=%d, nr=%d, on-surf GRF test at 1 pt: u err = %.3g\n',N,o.nr,u-uex)
        SAs(e) = sum(sn.w)+sum(w); %fprintf('area = %.12g\n',surfarea)
      end
      errs(e,c)=abs(u-uex);
    end        % ------------------------------------------
  end        % ===================
  
  % main convergence plot...
  dxs = sqrt(surfarea./Ns);   % defn of dx_typ for surface
  figure; %loglog(dxs,abs(SAs-surfarea),'o-'); hold on;
  loglog(dxs,errs(:,1),'.-','markersize',10); axis tight; hold on;
  loglog(dxs,errs(:,2:end),'+-');
  if fig==1, lab='(b)'; sc=0.1; else, lab='(c)'; sc=100; end  % label, offset
  loglog(dxs,sc*dxs.^o.p,'k--');     % plain order p (not pred 2p)
  %xlabel('$\Delta x_{typ}$','interpreter','latex');
  xlabel('$h$','interpreter','latex');  %'$\Delta x$'
  ylabel('pointwise error','interpreter','latex');
  h=legend('exterior',sprintf('on-surf. $n_r=%d$',nrs(2)),sprintf('on-surf. $n_r=%d$',nrs(3)),sprintf('on-surf. $n_r=%d$',nrs(4)),sprintf('order $%d$',o.p));
  set(h,'location','southeast','interpreter','latex');
  title(sprintf('%s surface quadrature convergence, p=%d',lab,o.p));
  %hold on; loglog(dxs,100*dxs.^(2*o.p),'m-'); % predicted 2p rate
  set(gcf,'paperposition',[0 0 3.5 3.5]);    % somehow black border in tex :(
  print('-depsc2',sprintf('spaceconv_p%d.eps',o.p));
end   %%%%%%%%%%%%%

if 1    % plot the ret tau used, making high-resolution smooth nodes for plot...
  so.np = 15; so.mp = 10; so.orig=0; o.p = 8;   % new params
  [s N] = create_panels('torus',so,o);  % copied from above
  [x nx w] = getallnodes(s);  surfarea = sum(w);  % use large np,mp to get S.A.
  t.x = s{1}.x(:,1);  % on-surf, as above
  tret = ttarg - dists(t.x,x);
  fprintf('min dist of src for data (xs) from surface: %.3g\n',min(dists(xs,x)))
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);  % eval dens, copied from above
  retsig = -fn; rettau = f; rettaut = ft;     % (note retarded dens dep on targ)
  showsurffunc(s,rettau); colorbar off; hold on;
  view(35,5); %view(60,40);
  plot3(t.x(1),t.x(2),t.x(3),'k.','markersize',10);
  plot3(xext(1),xext(2),xext(3),'k.','markersize',10);
  plot3(xs(1),xs(2),xs(3),'b.','markersize',10);
  axis off; text(-1,-1,1,'(a)');
  set(gcf,'paperposition',[0 0 4 3]);
  print -dpng -r600 spaceconv_rettau.png
  system('convert -trim spaceconv_rettau.png eps2:spaceconv_rettau.eps')
end


if 0 & side~=0   % (old plot for off-surf case only)
  fprintf('min dist of trg (t.x) from surface: %.3g\n',min(dists(t.x,x)))
  showsurffunc(s,rettau); title('GRF: retarded tau');
  view(40,30);
  hold on; plot3(t.x(1),t.x(2),t.x(3),'k.','markersize',10);
  %plot3(xext(1),xext(2),xext(3),'kx','markersize',20);
  colorbar off;
end

%%%%%%%%%%%%%%%%%% Indep warm-up tests...
if 0  % check expected composite gauss quad order is 2p, on [0,1]. Yes.
  w0 = 5.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t);
  Iex = T(1)-T(0);   % let's integrate Tt
  p = 5;
  [x w] = gauss(p);
  ns = 2:2:10; errs = nan*ns;  % n=# panels in [0,1]
  for i=1:numel(ns), n=ns(i)
    I = 0;
    for j=1:n, t=((x-1)/2+j)/n; I = I+(1/2/n)*w*Tt(t); end
    errs(i) = abs(I-Iex);
  end
  figure; loglog(ns,errs,'+-'); hold on; loglog(ns,2e-4*ns.^(-2*p),'r-');
  axis tight; xlabel('# panels'); ylabel('error');
end
%%%%%%%%%%%%%%%%%%%%%%