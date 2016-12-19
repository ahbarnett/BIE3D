% test t-domain GRF for wave equation outside/inside/on a torus surface.
% With Hagstrom+Greengard.
% Barnett 12/15/16

clear
so.a=1; so.b = 0.5; [s N] = create_panels('torus',so,[]); % surf: default # pans

% surf data...
%t0 = 1.0; T = @(t) exp(-(t/t0).^2/2); Tt = @(t) (-t/t0^2).*T(t); % data t-func
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t);
eps = 1e-5; t = 0.3;       % test time
fprintf('error in time deriv: %.3g\n',(T(t+eps)-T(t-eps))/(2*eps) - Tt(t))
clear t
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside

side = 0;  % -1,0,1: choose nature of GRF test pt (a fake 1-pt panel)

if side==1, t.x = [1.3;0.1;0.8];          % exterior...
elseif side==-1, t.x = [-.8;.1;.2];       % or interior
else, k=57;j=1; t.x = s{k}.x(:,j); end    % or on-surf
t.N = 1; t.t = 2.1;                       % test target time
[x nx] = getallnodes(s);
r = sqrt(sum(bsxfun(@minus,t.x,x).^2,1));   % dists from test pt to surf pts
tret = t.t - r;   % retarded times on surf rel to GRF test pt

if side~=0          % ------------------- off-surf 
  % now eval retarded sig, tau, tau'...
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  %showsurffunc(s,retsig); title('GRF: retarded sigma'); showsurffunc(s,rettau); title('GRF: retarded tau'); showsurffunc(s,rettaut); title('GRF: retarded tau_t');

  % do u(x,t) = S.sigma + (wave eqn D).tau :
  u = LapSeval_panels(t,s,retsig) + LapDeval_panels(t,s,rettau) + RetDeval_panels(t,s,rettaut)
  uex = data_ptsrc(xs,T,Tt,t.t,t.x) * (side==1)           % zero inside, f out
  fprintf('N=%d, off-surf GRF test at 1 pt: u err = %.3g\n', N, u-uex)
  % error 1e-8

else               % ------------- on-surf GRF w/ local aux nodes, pan k, pt j
  o.nr = 8; o.nt = 2*o.nr; s = add_panels_auxquad(s,o);    % setup aux quad
  s = add_panels_timequad(s);                              % local aux t-delays
  kill = [k; s{k}.nei];     % pan inds to remove from far list
  sf = s; sf(kill) = [];    % far pans - need fix internal inds in sf ???
  [x nx] = getallnodes(sf);
  r = sqrt(sum(bsxfun(@minus,t.x,x).^2,1));   % dists from test pt to surf pts
  tret = t.t - r;   % retarded times on surf rel to GRF test pt
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);      % far
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  u = LapSeval_panels(t,sf,retsig) + LapDeval_panels(t,sf,rettau) + RetDeval_panels(t,sf,rettaut);
  % make fake near panel w/ all aux nodes and hack the time delays for now...
  sn.x = s{k}.auxnodes(:,:,j);  sn.nx = s{k}.auxnormals(:,:,j);
  sn.w = s{k}.auxwei(:,j)';           % must be col vec
  sn.N = numel(sn.w);       % finish the fake aux src pan
  tret = t.t - s{k}.auxdelays(:,j)';   % near panel's aux eval times, col vec
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,sn.x,sn.nx);      % near
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  u = u + LapSeval_panels(t,sn,retsig) + LapDeval_panels(t,sn,rettau) + RetDeval_panels(t,sn,rettaut);
  uex = data_ptsrc(xs,T,Tt,t.t,t.x) / 2;       % on-surf GRF
  fprintf('N=%d, on-surf GRF test at 1 pt: u err = %.3g\n', N, u-uex)
end

% ** note conflict field btw s.t - time for nodes, vs coord chart locs in R2!
