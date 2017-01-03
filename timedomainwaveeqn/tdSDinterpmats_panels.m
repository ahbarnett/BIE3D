function [Sret Dret] = tdSDinterpmats_panels(tpan,span,Linfo,intinfo)
%
% [Sret Dret] = tdSDinterpmats_panels(tpan,span,Linfo,intinfo)
%
% eval p^2-by-NT matrices which apply retarded S,D ops from dens hist grids
%
% where N is # dofs in all src pans.
%
% tpan - target "panel" (ie has fields t.x, t.N)
%       If is same as a source panel in span, tpan must have auxnodes.
% Linfo - interpolation struct for the one std panel (from setup_auxinterp).
% intinfo - time-interp struct, has fields:
%         n = # history steps, dt = timestep, m = interp order.
%
% tpan single panel for now. Includes aux node close & self-eval.
%
% See also: TEST_TDGRF_INTERP.m which tests this (near bottom).

% Barnett 12/29/16
if nargin==0, test_tdSDinterpmats_panels; return; end

M = size(getallnodes({tpan}),2);          % # targets
t = tpan;
if numel(span)==1, span = {span}; end
% ** idea: loop through all span, tra
% If it's a self-or-nei of t, then treat each targ node in t separately,
% since they have own aux nodes (to get spatial sing correct).
% Will need to make dummy targ panel w/ 1 node, loop over t's targ nodes.
n = intinfo.n;
Sret = []; Dret = [];
for q=1:numel(span), s = span{q};    % loop over src pans in right order
  r = relatedpanel(t,s);
  if r==0
    [Sq Dq] = tdSDinterpmats_panelpair(t,s,intinfo);
  else    % s is self or nei of t
    nN = s.N*n;
    Sq = nan(t.N,nN); Dq = Sq;
    for j=1:t.N             % loop over targs and write each as row of Sq
      tj.x = t.x(:,j); tj.N = 1;    % this targ pt as own struct
      i = Linfo.auxindsbytarg{r}{j};  % indices in full list of aux nodes
      saux.x = t.auxnodes(:,i);   % here i indexes the last 2 dims naux*N
      saux.nx = t.auxnormals(:,i);
      saux.w = t.auxwei(i);
      saux.N = numel(saux.w);
      [Sa Da] = tdSDinterpmats_panelpair(tj,saux,intinfo); % use saux as src pan
      Ltjsaux = Linfo.Lbytarg{r}{j};     % see: setup_auxinterp
      Sq(j,:) = reshape(reshape(Sa,[n saux.N]) * Ltjsaux, [1 nN]);
      Dq(j,:) = reshape(reshape(Da,[n saux.N]) * Ltjsaux, [1 nN]);
    end
  end
  Sret = [Sret, Sq]; Dret = [Dret, Dq];   % stack each pan as block col
end

%%%%%%
function [Sret Dret] = tdSDinterpmats_panelpair(t,s,o)
% Inputs:
% t - target panel struct with: t.x - 3xM target locs
% s - source panel struct with: s.x, s.nx - 3xN locs and normal, s.w 1xN weights
% o - interpolation info struct with fields: n (# history steps), dt, m.
% Outputs:
% Sret,Dret - M-by-Nn sparse matrices, each row of which is a vector to apply
%             appropriately retarded SLP or DLP to density history vectors
%             (ordered with time fast, nodes slow) for that row's target.
%             Thus these matrices can do true matvecs against dens history.
% Barnett 12/29/16
n = o.n;
M = size(t.x,2); N = numel(s.w);        % # targs, # srcs
[S D Dp] = tdSDmats(t.x,s.x,s.nx,s.w);  % spatial quadr mats, each is MxN
delays = dists(s.x,t.x);                % pt pairwise time delays >0, transpose
[jmax jmin A Ap] = interpmat(-delays(:),o.dt,o.m);  % Tom's weights (1 row per delay, ordered fast over sources, slow over targs)
joff = jmin+o.n-1;         % padding on the ancient side
if joff<0, error('interp requesting too ancient history!'); end
ii = []; jj = []; aa = []; dd = [];  % to build the sparse mats
for k=1:M     % loop over targs
  [j i a] = find(A((1:N)+(k-1)*N,:)');    % j is time inds, i is src inds (slow)
  aa = [aa; a.*S(k,i)'];                  % SLP spatial kernel & quadr wei
  dd = [dd; a.*D(k,i)'];                  % value part of DLP
  ii = [ii; k*ones(size(a))];
  jj = [jj; joff+j+n*(i-1)];              % time indices in the Nn vector
end
Sret = sparse(ii,jj,aa,M,N*n);
Dret = sparse(ii,jj,dd,M,N*n);            % (only the value part so far)
% *** todo: find neater way to build this without repeating the find()...
% (issue is want sparsity pattern that includes A and Ap...)
ii = []; jj = []; dd = [];  % to build the deriv part of DLP
for k=1:M     % loop over targs
  [j i a] = find(Ap((1:N)+(k-1)*N,:)');   % j is time inds, i is src inds (slow)
  dd = [dd; a.*Dp(k,i)'];                 % deriv part of DLP
  ii = [ii; k*ones(size(a))];
  jj = [jj; joff+j+n*(i-1)];              % time indices in the Nn vector
end
Dret = Dret + sparse(ii,jj,dd,M,N*n);  % may have slightly different patterns


%%%%%
function test_tdSDinterpmats_panels   % do off/on surf wave eqn GRF test,
% taken from test_tdGRF_interp.  Barnett 1/1/17
side = 0;    %  GRF test:   1 ext, 0 on-surf
dt = 0.1;   % timestep
m = 4;      % control time interp order (order actually m+2)

so.a=1; so.b=0.5; o.p=6;
[s N] = create_panels('torus',so,o); % surf: default # pans
[x nx w] = getallnodes(s);
distmax = 4.0;       % largest dist from anything to anything
n = ceil(distmax/dt);

if side==1
  t.N = 1; t.x = [1.3;0.1;0.8];    % single test targ pt, exterior...
  Linfo = [];              % spatial interp info
else
  o.nr = 8; o.nt = 2*o.nr;     % first add aux quad to panels: aux quad orders
  s = add_panels_auxquad(s,o);
  t = s{57};          % on-surf targ is a whole panel
  Linfo = setup_auxinterp(s{1}.t,o);  % std spatial interp to aux quad
end
ttarg = 0.0;          % test target time (avoids "t" panel field conflict)

% surf data for GRF...
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t); % data src t-func, tested
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside
% eval sig, tau on {n history grid} x {N bdry nodes}
tt = dt*(-n+1:0); ttt = repmat(tt,[1 N]);
xx = kron(x,ones(1,n)); nxx = kron(nx,ones(1,n));   % ttt,xx,nxx spacetime list
[f,fn] = data_ptsrc(xs,T,Tt,ttt,xx,nxx);       % output ft unused
sighist = -fn; tauhist = f;  % col vecs, ext wave eqn GRF: u = D.u - S.un

[Starg,Dtarg] = tdSDinterpmats_panels(t,s,Linfo,struct('n',n,'dt',dt,'m',m));

u = Starg*sighist + Dtarg*tauhist;
uex = data_ptsrc(xs,T,Tt,ttarg,t.x);     % what ext GRF should give
if side==0, uex=uex/2; end   % on-surf principal value
fprintf('N=%d, dens-interp ext GRF test: u err = %.3g\n', N, max(abs(u-uex)))
