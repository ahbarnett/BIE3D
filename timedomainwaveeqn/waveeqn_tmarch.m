% First baby timesteps in t-domain wave eqn BIE: known exterior Dirichlet BVP.
% Barnett 1/4/17 w/ Hagstrom, Greengard.

clear
dt = 0.1;   % timestep
m = 4;      % control time interp order (order actually m+2)

so.a=1; so.b=0.5; o.p=6;
[s N] = create_panels('torus',so,o); % surf: default # pans
[x nx w] = getallnodes(s);
distmax = 4.0;       % largest dist from anything to anything
n = ceil(distmax/dt);

SDload = 1;
if SDload
  load SDtarg_torus_p6_m4_dt01      % precomputed 2GB (took 80 sec)
else
  o.nr = 8; o.nt = 2*o.nr;     % first add aux quad to panels: aux quad orders
  s = add_panels_auxquad(s,o);
  Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
  [Starg,Dtarg] = tdSDinterpmats_panels(t,s,Linfo,struct('n',n,'dt',dt,'m',m));
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
tret = ttarg - dists(t.x,x);     % retarded times of surf nodes rel to test pt
[jmax,jmin,a,ap] = interpmat(tret,dt,m);    % Tom's coeffs (1 row per tret)
joff = jmin+n-1;         % padding on the ancient side
if joff<0, error('u_test eval vec requesting too ancient history!'); end
a = [zeros(N,joff), a, zeros(N,-jmax)];  % pad to width n, preserve sparsity
ap = [zeros(N,joff), ap, zeros(N,-jmax)];
[S D Dp] = tdSDmats(t.x,x,nx,w);  % each is 1xN
Stest = a'.*repmat(S,[n 1]);      % coeff vectors packed as nxN matrices
Dtest = a'.*repmat(D,[n 1]) + ap'.*repmat(Dp,[n 1]);
Stest = Stest(:)'; Dtest = Dtest(:)';  % the u-test-eval row vecs

utest = Stest*sighist + Dtest*tauhist;
uex = data_ptsrc(xs,T,Tt,0.0,t.x);     % exact u at test pt
fprintf('test that the u eval vectors (SDtest) work: %.3g\n', utest-uex)

% pull out current-time matrix (last time index not first!):
Snow = Starg(:,n:n:end); Dnow = Dtarg(:,n:n:end);  % NxN
% rep is u = D.mu:
[Lnow Unow pnow] = lu(speye(N)/2 + Dnow,'vector');   % direct factor on-surf sys
%rhs = randn(N,1); mu = Unow\(Lnow\rhs(pnow)); norm(Dnow*mu + mu/2 - rhs) % check direct solve via LU
% [set up sparse to do extrap on each node separately (only for pred-corr)]

% Dirichlet data for BVP: t-dep pulse func for interior pt src...
t0=6; s0=1.0; T = @(t) exp(-0.5*(t-t0).^2/s0^2); Tt = @(t) -((t-t0)/s0^2).*T(t);
% (t0/s0 = 6 gives 1e-8 of start-up error if no time delay from src to surf)
%eps = 1e-5; tt = 4.3; fprintf('check Tt vs T: %.3g\n',(T(tt+eps)-T(tt-eps))/(2*eps) - Tt(tt)), clear tt

Ttot = 14.0;     % total time to evolve
jtot = ceil(Ttot/dt); 
verb = 1;
muhist = zeros(n*N,1);  % init dens hist
gs=nan(jtot,1); rs=gs; ms=gs; es=gs;   % to save sizes of things for later
if verb>1, figure; end
for j=1:jtot           % .... marching loop
  tj = j*dt;
  muhist = reshape(muhist,[n N]); muhist = muhist([2:n,n],:); % shuffle back 1
  muhist(end,:)=0; muhist = reshape(muhist,[n*N 1]); % ensure munow=0
  gnow = data_ptsrc(xs,T,Tt,tj,x,nx);                % Dirichlet data from src
  rhs = gnow - Dtarg*muhist;                         % RHS for "now" lin sys
  muhist(n:n:end) = Unow\(Lnow\rhs(pnow));           % write into "now" entries
  u = Dtest*muhist;                                  % eval at test pt
  uex = data_ptsrc(xs,T,Tt,tj,t.x);                  % known BVP soln at test pt
  fprintf('j=%d (tj=%.3g):\t err=%.3g \tu=%.6g \tuex=%.6g\n',j,tj,u-uex,u,uex)
  es(j)=u-uex; gs(j)=norm(gnow); rs(j)=norm(rhs); ms(j)=norm(muhist(n:n:end));
  %imagesc(reshape(muhist,[n N])); caxis([0 3]); colorbar; drawnow
  if verb>1, plot([gnow, rhs, muhist(n:n:end)], '.-'); title(sprintf('t=%.3g',tj)); legend('g','rhs','\mu_{now}'); drawnow, end
end                     % ....
if verb, figure; semilogy([gs rs ms abs(es)],'.-');
  legend('||g||','||rhs||','||\mu_{now}||','u err'); xlabel('timesteps');
  axis tight;
  title(sprintf('torus, non-osc pulse, \\delta t = %.3g',dt));
end





