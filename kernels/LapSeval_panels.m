function u = LapSeval_panels(tpan,span,dens,paninfo)
% LAPSEVAL_PANELS  evaluate Laplace SLP from surface in R3, incl on-surface
%
% u = LapSeval_panels(t,s,dens,paninfo)
%
% Called with no arguments does a self-test (see bottom of code for usage)
%
% s,t are src and targ panel cell arrays. t doesn't need valid t.inds (unused)
%  since the output dof ordering is built in order from the cell array t.
%  paninfo (only needed if self-eval detected) is struct containing
%  interpolation matrices L, etc; see setup_auxinterp.m.
%
%  does panel-by panel tests of whether self-eval; if true uses 3x3 patch
%  source auxiliary quadrature nodes, and does on-surface limit of value.
%
%  Note: can't handle on-surface targets that are not a complete smooth
%   quadr for a target panel; that would need recomputation of interp matrices!
%
% See also: TEST_LAPGRF for SLP test

% IDEAS: If this is called multiple times, prestore the aux kernel weights
% w/ this kernel, ie the local matrices.
% Make matrix version, which has dgemm (BLAS3) for ker mat * interp mat.

% Barnett 7/18/16. Much of code shared with LapDeval_panels

if nargin==0, test_LapSeval_panels; return; end
if nargin<4, paninfo = []; end

if numel(tpan)==1, tpan = {tpan}; end  % ensure cell arrays
if numel(span)==1, span = {span}; end
M = size(getallnodes(tpan),2);          % # targets
Np = numel(span);      % # src pans
n = span{1}.N;        % # nodes per src pan (assume same)
u = nan(M,1);
indoff = 0;        % track index in output array
for p = 1:numel(tpan), t = tpan{p}; % for each targ panel...
  tinds = indoff+(1:t.N);         % indices in output list
  spanselfk = find(samepanel(t,span));     % any self src pans? O(N), yuk
  if isempty(spanselfk)         % no self (use smooth for all src panels)
    [y ny w] = getallnodes(span);
    u(tinds) = LapSLPeval(t.x,y,w,dens);          % note ny absent!
  else                          % one src panel is self (same as targ)
    ss = span{spanselfk};       % the self src pan
    far = true(Np,1); far(ss.nei) = false; far(spanselfk) = false;
    [y ny w] = getallnodes({span{far}});    % note makes back a cell array
    farinds = logical(kron(far, ones(n,1)));  % assumes same n (nodes) for all src pans
    u(tinds) = LapSLPeval(t.x,y,w,dens(farinds));
    neardens = dens(~farinds);
    auxdens = zeros(1,ss.naux*t.N);  % step 1: interp dens to all aux nodes:
    for i=1:9                   % loop over self+nei 3x3 local (near) pans
      if i==1, k=spanselfk; else, k=span{spanselfk}.nei(i-1); end % pan index
      densk = dens(span{k}.inds);   % get density data for this pan
      auxdens(paninfo.auxinds{i}) = paninfo.L{i} * densk;
    end
    auxdens = reshape(auxdens,[ss.naux, t.N]);  % make easier to index
    unear = nan(t.N,1);            % step 2: apply kernel at aux nodes
    for i=1:t.N             % loop over targs; each has own set of aux nodes
      unear(i) = LapSLPeval(t.x(:,i),squeeze(ss.auxnodes(:,:,i)),ss.auxwei(:,i)',auxdens(:,i));
    end
    u(tinds) = u(tinds) + unear;
  end
  indoff = indoff + t.N;
end

function u = LapSLPeval(x,y,w,dens)
% Dumb dense Laplace SLP evaluator from sources w/ given weights, to targets
% Inputs:
% x targ (3-by-M), y src (3-by-N), w src weights (1-by-N)
%
% NB: fills M*N matrix, so assumes small cases
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
ir = 1./sqrt(d1.^2+d2.^2+d3.^2);   % 1/r matrix, size M*N
u = (1/4/pi) * (ir * (w(:).*dens));

%%%%%%%%%%%%%%%%%%%%%%
function test_LapSeval_panels    % tests just on-surface convergence of SLP
% ... not correctness of the formula; for that see test_LapGRF.m
so.a=1; so.b = 0.5; [s N] = create_panels('torus',so,[]);  % use default # pans
% non-const tau auxquad convergence (tests auxquad interp is good):
x = getallnodes(s);
dens = sin(1 + x(1,:) + 1.7*x(2,:) - 0.4*x(3,:))';  % smooth surf func, col vec
k = 57;       % which pan to test on-surf u vals
nrs = [4 6 8 10];
for i=1:numel(nrs), o.nr = nrs(i); o.nt = 2*o.nr;
  s = add_panels_auxquad(s,o);                      % overwrites
  fprintf('setup interp... '); Linfo = setup_auxinterp(s{1}.t,o);  % slow
  u = LapSeval_panels(s{k},s,dens,Linfo);  % potential should be -1/2 at all pts
  fprintf('SLP auxquad conv:   nr=%d\t u(1) = %.16g\n',o.nr, u(1)) % should conv
end
% gets about 1e-9 at nr=8
