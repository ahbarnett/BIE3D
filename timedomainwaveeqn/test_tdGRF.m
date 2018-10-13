% test t-domain GRF for wave equation outside/inside/on a torus surface.
% With Hagstrom+Greengard.
% Barnett 12/15/16

clear
wobbly = 1;
so.a = 1; b = 0.5;   % major radius and poloidal radius
o.p = 6;   % panel order
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
so.np = 12; so.mp = 8;  % default
%so.np = 18; so.mp = 12;
%so.np = 9; so.mp = 6;
[s N] = create_panels('torus',so,o); % surf: default # pans

% surf data...
%t0 = 1.0; T = @(t) exp(-(t/t0).^2/2); Tt = @(t) (-t/t0^2).*T(t); % data t-func
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t);
Tfunctest = 0; if Tfunctest
  eps = 1e-5; t = 0.3;       % test time
  fprintf('error in time deriv: %.3g\n',(T(t+eps)-T(t-eps))/(2*eps) - Tt(t))
  clear t
end
xs = [0.9;-0.2;0.1];   % src pt for data, must be well inside

side = 0;  % -1,0,1: choose nature of GRF test pt (a fake 1-pt panel)

if side==1, t.x = [1.3;0.1;0.8];          % exterior...
elseif side==-1, t.x = [-.8;-.3;.2];       % or interior
else, k=17;j=1; t.x = s{k}.x(:,j); end    % or on-surf
t.N = 1; ttarg = 2.1;              % test target time (avoid s.t conflict)
[x nx w] = getallnodes(s);
tret = ttarg - dists(t.x,x);   % retarded times on surf rel to GRF test pt

if side~=0          % ------------------- off-surf 
  % now eval retarded sig, tau, tau'...
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  %showsurffunc(s,retsig); hold on; plot3(t.x(1),t.x(2),t.x(3),'k.','markersize',20); title('GRF: retarded sigma');
  %showsurffunc(s,rettau); title('GRF: retarded tau'); showsurffunc(s,rettaut); title('GRF: retarded tau_t');
  % do u(x,t) = S.sigma + (wave eqn D).tau :
  u = LapSeval_panels(t,s,retsig) + LapDeval_panels(t,s,rettau) + RetDeval_panels(t,s,rettaut)
  [S D Dp] = tdSDmats(t.x,x,nx,w);
  v = S*retsig + D*rettau + Dp*rettaut;
  fprintf('max diff btw u eval @ targ via tdSDmats vs direct eval: %.3g\n',max(abs(v-u)))
  uex = data_ptsrc(xs,T,Tt,ttarg,t.x) * (side==1)         % zero inside, f out
  fprintf('N=%d, off-surf GRF test at 1 pt: u err = %.3g\n', N, u-uex)
  % error 4e-9 for wobbly=0 (2e-8 wobbly=1), at p=6 (np=12,mp=8), r=8

else               % ------------- on-surf GRF w/ local aux nodes, pan k, pt j
  o.nr = 8; o.nt = 2*o.nr; s = add_panels_auxquad(s,o);    % setup aux quad
  s = add_panels_timequad(s);                              % local aux t-delays
  kill = [k; s{k}.nei];     % pan inds to remove from far list
  sf = s; sf(kill) = [];    % far pans - need fix internal inds in sf ???
  [x nx] = getallnodes(sf);
  tret = ttarg - dists(t.x,x);   % retarded times on surf rel to GRF test pt
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);      % far
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  u = LapSeval_panels(t,sf,retsig) + LapDeval_panels(t,sf,rettau) + RetDeval_panels(t,sf,rettaut);
  % make fake near panel w/ all aux nodes and hack the time delays for now...
  sn.x = s{k}.auxnodes(:,:,j);  sn.nx = s{k}.auxnormals(:,:,j);
  sn.w = s{k}.auxwei(:,j)';           % must be col vec
  sn.N = numel(sn.w);       % finish the fake aux src pan
  tret = ttarg - s{k}.auxdelays(:,j)';   % near panel's aux eval times, col vec
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,sn.x,sn.nx);      % near
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  u = u + LapSeval_panels(t,sn,retsig) + LapDeval_panels(t,sn,rettau) + RetDeval_panels(t,sn,rettaut);
  uex = data_ptsrc(xs,T,Tt,ttarg,t.x) / 2;       % on-surf GRF
  fprintf('N=%d, on-surf GRF test at 1 pt: u err = %.3g\n', N, u-uex)
end
