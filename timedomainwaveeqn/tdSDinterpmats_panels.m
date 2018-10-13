function [Sret Dret Sdotret] = tdSDinterpmats_panels(tpan,span,Linfo,intinfo)
%
% [Sret Dret Sdotret] = tdSDinterpmats_panels(tpan,span,Linfo,intinfo)
%
%  eval p^2-by-NT matrices which apply retarded S,D, and S(d/dt) ops from dens
%  hist grids
%
%  where N is # dofs in all src pans.
%
% tpan - target "panel" (ie has fields t.x, t.N)
%       If is same as a source panel in span, tpan must have auxnodes.
% Linfo - interpolation struct for the one std panel (from setup_auxinterp).
% intinfo - time-interp struct, has fields:
%         n = # history steps, dt = timestep, m = interp order.
%
% tpan single panel for now. Includes aux node close & self-eval.
%
% Without arguments, does self-test.
%
% See also: TEST_TDGRF_INTERP.m which is an older driver (near bottom).

% Barnett 12/29/16 - 1/4/17. Sdotret 1/7/17. parfor 10/11/18
if nargin==0, test_tdSDinterpmats_panels; return; end

if numel(tpan)==1, tpan = {tpan}; end
M = size(getallnodes(tpan),2);          % # target nodes
if numel(span)==1, span = {span}; end
N = size(getallnodes(span),2);          % # source nodes
n = intinfo.n;
Si = cell(numel(tpan),1); Di = Si; Sdi = Si;  % cell array of structs, for parfor
parfor i=1:numel(tpan), t = tpan{i};   % ----- outer targ panel loop
  fprintf('targ pan #%d...\n',i)       % (tweaked so each i is indep task)
  coff = 0;                           % track column offset for sparse mat out
  nnzmax = ceil(20*N);                % allocation size for Si,Di targ-pan NNZ
  Si{i}.i = nan(nnzmax,1); Si{i}.j = Si{i}.i; Si{i}.v = Si{i}.i; Si{i}.ptr = 0;
  Di{i} = Si{i}; Sdi{i} = Si{i};
  for q=1:numel(span), s = span{q};   % loop over src pans in right order
    r = relatedpanel(t,s);
    if r==0    % s is unrelated, ie far from t
      [Sq Dq Sdq] = tdSDinterpmats_panelpair(t,s,intinfo);
    else    % s is self or nei of t
      nN = s.N*n;
      Sq = nan(t.N,nN); Dq = Sq; Sdq = Sq;
      for j=1:t.N             % loop over targs and write each as row of Sq
        tj = struct();
        tj.x = t.x(:,j); tj.N = 1;    % this targ pt as own struct
        k = Linfo.auxindsbytarg{r}{j};  % indices in full list of aux nodes
        saux = struct();
        saux.x = t.auxnodes(:,k);   % here i indexes the last 2 dims naux*N
        saux.nx = t.auxnormals(:,k);
        saux.w = t.auxwei(k);
        saux.N = numel(saux.w);         % we now use saux as src pan...
        [Sa Da Sda] = tdSDinterpmats_panelpair(tj,saux,intinfo);
        Ltjsaux = Linfo.Lbytarg{r}{j};     % see: setup_auxinterp
        Sq(j,:) = reshape(reshape(Sa,[n saux.N]) * Ltjsaux, [1 nN]);
        Dq(j,:) = reshape(reshape(Da,[n saux.N]) * Ltjsaux, [1 nN]);
        Sdq(j,:) = reshape(reshape(Sda,[n saux.N]) * Ltjsaux, [1 nN]);
      end
    end
    %  dump each src blk into sparse lists for this i'th targ panel...
    [ii jj vv] = find(Sq); nh=numel(ii); hh = Si{i}.ptr+(1:nh); Si{i}.ptr = Si{i}.ptr+nh;
    Si{i}.i(hh) = ii;  Si{i}.j(hh) = jj+coff;  Si{i}.v(hh) = vv;
    [ii jj vv] = find(Dq); nh=numel(ii); hh = Di{i}.ptr+(1:nh); Di{i}.ptr = Di{i}.ptr+nh;
    Di{i}.i(hh) = ii;  Di{i}.j(hh) = jj+coff;  Di{i}.v(hh) = vv;
    [ii jj vv] = find(Sdq); nh=numel(ii); hh = Sdi{i}.ptr+(1:nh); Sdi{i}.ptr = Sdi{i}.ptr+nh;
    Sdi{i}.i(hh) = ii;  Sdi{i}.j(hh) = jj+coff;  Sdi{i}.v(hh) = vv;
    coff = coff + size(Sq,2);
  end
end                                % ------

t0=tic;
S.i = []; S.j = []; S.v = []; S.ptr = 0; D = S; Sd = S;  % full sparse lists
roff = 0;                           % track row offset in sparse mat out
for i=1:numel(tpan)  % copy each targ pan into full sparse lists...(use parfor?)
  % dump each targ blk row into sparse lists... (don't append; too slow!)
  hh = S.ptr+(1:Si{i}.ptr); S.ptr=S.ptr+Si{i}.ptr;  % inds in final sparse list
  S.i(hh)=Si{i}.i+roff; S.j(hh)=Si{i}.j; S.v(hh)=Si{i}.v;
  hh = D.ptr+(1:Di{i}.ptr); D.ptr=D.ptr+Di{i}.ptr;  % "
  D.i(hh)=Di{i}.i+roff; D.j(hh)=Di{i}.j; D.v(hh)=Di{i}.v;
  hh = Sd.ptr+(1:Sdi{i}.ptr); Sd.ptr=Sd.ptr+Sdi{i}.ptr; % "
  Sd.i(hh)=Sdi{i}.i+roff; Sd.j(hh)=Sdi{i}.j; Sd.v(hh)=Sdi{i}.v;
  roff = roff + tpan{i}.N;
end
fprintf('stack cells into single sparse lists: %.3g s\n',toc(t0))

if exist('fsparse')~=3  % unlike panelpair case, matlab assemble can be beaten:
  t0=tic;
  Sret = sparse(S.i,S.j,S.v,M,n*N);     % build each matrix in one go
  Dret = sparse(D.i,D.j,D.v,M,n*N);
  Sdotret = sparse(Sd.i,Sd.j,Sd.v,M,n*N);
  fprintf('matlab sparse build: %.3g s\n',toc(t0))
else
  t0=tic;   % speed up w/ multithreaded stenglib/Fast/fsparse.c MEX :
  % it's possible 'nosort' speeds up build or spmatvec when added here... (nope)
  Sret = fsparse(S.i,S.j,S.v,[M,n*N,numel(S.i)]);  % build each matrix in one go
  Dret = fsparse(D.i,D.j,D.v,[M,n*N,numel(D.i)]);
  Sdotret = fsparse(Sd.i,Sd.j,Sd.v,[M,n*N,numel(Sd.i)]);
  fprintf('stenglib fsparse build: %.3g s\n',toc(t0))
end
  
% Notes: spreplace here was too slow:
% https://www.mathworks.com/matlabcentral/answers/69528-sparse-matrix-more-efficient-assignment-operation
% Also, repeated append too slow


%%%%%%
function [Sret Dret Sdotret] = tdSDinterpmats_panelpair(t,s,o)
% Inputs:
% t - target panel struct with: t.x - 3xM target locs
% s - source panel struct with: s.x, s.nx - 3xN locs and normal, s.w 1xN weights
% o - interpolation info struct with fields: n (# history steps), dt, m.
% Outputs:
% Sret,Dret,Sdotret - M-by-Nn sparse matrices, each row of which is a vector
%             to apply appropriately retarded SLP, DLP, or SLP(d/dt .) to
%             density history vectors (ordered with time fast, nodes slow) for
%             that row's target.
%             Thus these matrices can do true matvecs against dens history.
% Barnett 12/29/16. Sdot 1/6/17. Removed appending to lists & find 10/10/18
n = o.n;
M = size(t.x,2); N = numel(s.w);        % # targs, # srcs
[S D Dp] = tdSDmats(t.x,s.x,s.nx,s.w);  % spatial quadr mats, each is MxN
delays = dists(s.x,t.x);                % pt pairwise time delays >0, transpose
[~,jmin,A,Ap] = interpmat(-delays(:),o.dt,o.m);  % Tom's weights (1 col per delay, row ordering fast over sources, slow over targs)
joff = jmin+o.n-1;                      % padding on the ancient side
if joff<0, error('interp requesting too ancient history!'); end
% get j (time) & i (src-targ pair) inds of all nonzero els...
[jj ii] = find(abs(A')+abs(Ap')>1e-12); % joint sparsity. transp for ii contig
aa = A(sub2ind(size(A),ii,jj));         % get all interp mat el vals
aap = Ap(sub2ind(size(A),ii,jj));
dd = 0*aa;
ir = 0*ii; jr = ir;                     % ind lists for retarted M*Nn mats
brk = [0; find(diff(floor((ii-1)/N))); numel(ii)];  % ind breakpts btw targs
for k=1:M                               % loop over targs
  i = (brk(k)+1 : brk(k+1))';           % col of inds of sparse inds for kth targ
  %i = find(ii>(k-1)*N & ii<=k*N)       % (same math, was slower)
  iii = ii(i)-(k-1)*N;                  % indices of src nodes
  dd(i) = aa(i).*D(k,iii)' + aap(i).*Dp(k,iii)';  % value & deriv parts of DLP
  % (NB since we're done with aa and aap for DLP, can now change them!)
  aa(i) = aa(i).*S(k,iii)';             % SLP spatial kernel & quadr wei
  aap(i) = aap(i).*S(k,iii)';           % SLP spatial acting on d/dt dens
  ir(i) = k;                            % all same targ index for ret mats
  jr(i) = joff+jj(i)+n*(iii-1);         % time indices within the Nn vector
end
if 1 | exist('fsparse')~=3    % native assemble - is fastest, always use
  Sret = sparse(ir,jr,aa,M,N*n);
  Sdotret = sparse(ir,jr,aap,M,N*n);
  Dret = sparse(ir,jr,dd,M,N*n);
else % stenglib, is strangely 3x slower than matlab here, so don't use.
  Sret = fsparse(ir,jr,aa,[M,N*n,numel(ir)]);
  Sdotret = fsparse(ir,jr,aap,[M,N*n,numel(ir)]);
  Dret = fsparse(ir,jr,dd,[M,N*n,numel(ir)]);
end

%%%%%
function test_tdSDinterpmats_panels   % do off/on surf wave eqn GRF S,D test,
% taken from test_tdGRF_interp.  Barnett 1/1/17. Not doesn't test Sdot.
side = 0;    %  GRF test:   1 ext, 0 on-surf
bigtest = 0;   % use all pans as on-surf targs (takes 2 mins)
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
  if bigtest, t = s; else t = s{57}; end % on-surf targ, 1 or more pans s(57:86)
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
sighist = -fn; tauhist = f;    % col vecs, ext wave eqn GRF: u = D.u - S.un

t0=tic; %profile clear; profile on
[Starg,Dtarg] = tdSDinterpmats_panels(t,s,Linfo,struct('n',n,'dt',dt,'m',m));
fprintf('total build time %.3g s\n',toc(t0)), %profile off; profile viewer

%[ii jj] = find(Dtarg); numel(ii)/prod(size(Dtarg)) % check sparsity
% tic, [ii jj aa] = find(Dtarg); toc % timing test

tic; u = Starg*sighist + Dtarg*tauhist; fprintf('t-step matvec in %.3g s\n',toc)
xtarg = getallnodes(t);
uex = data_ptsrc(xs,T,Tt,ttarg,xtarg);     % what ext GRF should give
if side==0, uex=uex/2; end   % on-surf principal value
fprintf('N=%d, dens-interp ext GRF test: max u err = %.3g\n',N,max(abs(u-uex)))
%u, uex
%whos   % 30MB, 0.6 sec per pan -> 3GB, 1 min for whole surf.
%keyboard

% todo: devise some test of Sdottarg, eg seeing if close to Starg applied to
% independently t-deriv'ed mu.

