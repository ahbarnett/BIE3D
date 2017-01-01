% test t-domain GRF for wave equation outside/inside/on a torus surface,
% from uniform-time-grid density on panel nodes.
% This is basically test_tdGRF but w/ t-interp from sig,tau on uniform t grid.
% Limitations: Tom's coeff mats are basically dense since all dists bundled.
% side=1 (ext) only for now.
% Barnett 12/28/16 - 12/29/12 for Hagstrom+Greengard project.

clear; verb = 1;
so.a=1; so.b=0.5; o.p=6;
[s N] = create_panels('torus',so,o); % surf: default # pans
diam = 2*(so.a+so.b);                       % specific to torus shape

% surf data...
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t); % data src t-func, tested
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside
[x nx w] = getallnodes(s);
dt = 0.1;   % timestep
m = 4;      % control time interp order (order actually m+2)

t.N = 1; t.x = [1.3;0.1;0.8];    % single test targ pt, exterior...
ttarg = 0.0;          % test target time (avoids "t" panel field conflict)
tret = ttarg - dists(t.x,x);     % retarded times

n = ceil(10 + m + max([-tret, diam])/dt);   % # history t-steps to store
[jmax jmin a ap] = interpmat(tret,dt,m);    % Tom's coeffs (1 row per tret)
fprintf('jrange = [%d,%d] in n=%d time history steps\n',jmin,jmax,n)
N = size(x,2);                              % # src nodes
a = [zeros(N,jmin+n-1), a, zeros(N,-jmax)];  % pad to width n, preserve sparsity
ap = [zeros(N,jmin+n-1), ap, zeros(N,-jmax)];

% eval sig, tau on {n history grid} x {N bdry nodes}
tt = dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,T,Tt,ttt,xx,nxx);       % output ft unused
sighist = -reshape(fn,[n N]); tauhist = reshape(f,[n N]);  % ext GRF: D.u - S.un

if verb % check interp matches direct eval of densities ret by eval at targ...
  sigri = sum(a.*sighist',2);                  % retarded dens got via interp
  tauri = sum(a.*tauhist',2); tautri = sum(ap.*tauhist',2);
  [f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);
  retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
  fprintf('interp err, ret dens from dens hist : sig %.2g, tau %.2g, tau'' %.2g\n',max(abs(sigri-retsig)), max(abs(tauri-rettau)), max(abs(tautri-rettaut)))
end   % ...we notice loss of 2 digits when interp to give a t-deriv :(

% build vectors which eval u(t=0,x=targ) given dens history data:
[S D Dp] = tdSDmats(t.x,x,nx,w);  % each is 1xN
Starg = a'.*repmat(S,[n 1]);      % coeff vectors packed as nxN matrices
Dtarg = a'.*repmat(D,[n 1]) + ap'.*repmat(Dp,[n 1]);
Starg = Starg(:)'; Dtarg = Dtarg(:)';  % row vecs

% or test the routine to build them
t2.x = [t.x, t.x+[1;0;0]];    % 2-target test "panel", checks multi-targ case
[Starg,Dtarg] = tdSDinterpmats_panels(t2,s,[],struct('n',n,'dt',dt,'m',m));
Starg = Starg(1,:); Dtarg = Dtarg(1,:); % keep only 1st target

% test them on dens history data, via u(x,t) = S.sigma + (wave eqn D).tau :
u = dot(Starg,sighist(:)) + dot(Dtarg,tauhist(:));
uex = data_ptsrc(xs,T,Tt,ttarg,t.x);     % what ext GRF should give
fprintf('N=%d, dens-interp ext GRF test at 1 pt: u err = %.3g\n', N, u-uex)

% observe: error a couple digits better than raw dens interp error,
% even if t-interp is low-order and dt large. :)
