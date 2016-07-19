function u = LapDeval_panels(tpan,span,dens,paninfo)
% LAPDEVAL_PANELS  evaluate Laplace DLP from surface in R3, incl on-surface
%
% u = LapDeval_panels(t,s,dens,paninfo)
%
% Called with no arguments does a self-test (see bottom of code for usage)
%
% s,t are src and targ panel cell arrays. t doesn't need valid t.inds (unused)
%  since the output dof ordering is built in order from the cell array t.
%  paninfo (optional) is struct containing interpolation matrices L, etc; see
%  setup_auxinterp.m.
%  (If paninfo is empty for self-evaluation, a const dens=1 is used, for DLP
%   testing only!)
%
%  does panel-by panel tests of whether self-eval; if true uses 3x3 patch
%  source auxiliary quadrature nodes, and does on-surface limit (no jump
%  relation).
%
%  Note: can't handle on-surface targets that are not a complete smooth
%   quadr for a target panel; that would need recomputation of interp matrices!

% IDEAS: If this is called multiple times, prestore the aux kernel weights
% w/ this kernel, ie the local matrices.
% Make matrix version, which has dgemm (BLAS3) for ker mat * interp mat.

% Barnett 7/15/16 - 7/17/16
if nargin==0, test_LapDeval_panels; return; end
if nargin<4, paninfo = []; end

if numel(tpan)==1, tpan = {tpan}; end  % ensure cell arrays
if numel(span)==1, span = {span}; end
tp = horzcat(tpan{:}); M = sum(horzcat(tp.N));  % # targets
Np = numel(span);      % # src pans
n = span{1}.N;        % # nodes per src pan (assume same)
u = nan(M,1);
indoff = 0;        % track index in output array
for p = 1:numel(tpan), t = tpan{p}; % for each targ panel...
  tinds = indoff+(1:t.N);         % indices in output list
  spanselfk = find(samepanel(t,span));     % any self src pans? O(N), yuk
  if isempty(spanselfk)         % no self (use smooth for all src panels)
    p = horzcat(span{:});       % make struct array (slow copying crap?)
    w = horzcat(p.w); y = horzcat(p.x); ny = horzcat(p.nx);
    u(tinds) = LapDLPeval(t.x,y,ny,w,dens);
  else                          % one src panel is self (same as targ)
    ss = span{spanselfk};       % the self src pan
    far = true(Np,1); far(ss.nei) = false; far(spanselfk) = false;
    p = horzcat(span{far});     % slow?
    w = horzcat(p.w); y = horzcat(p.x); ny = horzcat(p.nx);
    farinds = logical(kron(far, ones(n,1)));  % assumes same n (nodes) for all src pans
    u(tinds) = LapDLPeval(t.x,y,ny,w,dens(farinds));
    neardens = dens(~farinds);
    if isempty(paninfo)   % TESTING CODE FOR tau=1 ONLY!
      unear = nan(t.N,1);
      for i=1:t.N             % loop over targs; each has own set of aux nodes
        auxdens = ones(ss.naux,1);  % fix tau=1 on this targ's aux nodes
        unear(i) = LapDLPeval(t.x(:,i),squeeze(ss.auxnodes(:,:,i)),squeeze(ss.auxnormals(:,:,i)),ss.auxwei(:,i)',auxdens);
      end
    else                  % the actual code
      auxdens = zeros(1,ss.naux*t.N);  % step 1: interp dens to all aux nodes:
      for i=1:9                   % loop over self+nei 3x3 local (near) pans
        if i==1, k=spanselfk; else, k=span{spanselfk}.nei(i-1); end % pan index
        densk = dens(span{k}.inds);   % get density data for this pan
        auxdens(paninfo.auxinds{i}) = paninfo.L{i} * densk;
       % if i==3  % ** DEBUGGING: check overlays torus plot panel 49(3)
       %   aa = reshape(ss.auxnodes,[3 ss.naux*t.N]);
       %   l = aa(:,paninfo.auxinds{i});  % aux nodes
       %   plot3(l(1,:),l(2,:),l(3,:),'y.');
       % end
      end
      auxdens = reshape(auxdens,[ss.naux, t.N]);  % make easier to index
      unear = nan(t.N,1);            % step 2: apply kernel at aux nodes
      for i=1:t.N             % loop over targs; each has own set of aux nodes
        unear(i) = LapDLPeval(t.x(:,i),squeeze(ss.auxnodes(:,:,i)),squeeze(ss.auxnormals(:,:,i)),ss.auxwei(:,i)',auxdens(:,i));
      end
    end
    u(tinds) = u(tinds) + unear;
  end
  indoff = indoff + t.N;
end

function u = LapDLPeval(x,y,ny,w,dens)
% Dumb dense Laplace DLP evaluator from sources w/ given weights, to targets
% Inputs:
% x targ (3-by-M), y src, ny src normals (both 3-by-N), w src weights (1-by-N)
%
% NB: fills M*N matrix, so assumes small cases
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
rr = d1.^2+d2.^2+d3.^2;   % M*N
ddotn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
A =  ddotn ./ (sqrt(rr).*rr);  % kernel matrix w/o 1/4pi or src quad weights
u = (1/4/pi) * (A * (w(:).*dens));

%%%%%%%%%%%%%%%%%%%%%%
function test_LapDeval_panels    % tests some off- and on-surface stuff
so.a=1; so.b = 0.5; [s N] = create_panels('torus',so,[]);  % use default # pans
dens = ones(N,1);
t.x = [0.9;-0.2;0.1]; t.N = 1; % targ pt (as fake panel), must be inside
u = LapDeval_panels(t,s,dens);  % potential should be -1
fprintf('interior DLP tau=1 err from -1 (tests far off-surf): %.3g\n',u+1)

o.nt = 16; o.nr = 8;   % Test on-surf singular quadr: # aux nodes in theta,rho
s = add_panels_auxquad(s,o);
tt = s{1}.t;         % any panel will do to pass in the smooth param grid
%Linfo = [];   % for tau=1 test bypassing interp
Linfo = setup_auxinterp(tt,o);  % full test (interpolates 1: not much of test!)
k = 57;               % which pan to use as on-surf target pan
u = LapDeval_panels(s{k},s,dens,Linfo);  % potential should be -1/2 at all pts
fprintf('on-surface DLP tau=1 err from -1/2 (max over 1 pan): %.3g\n',max(abs(u+0.5)))

% non-const tau auxquad convergence (tests auxquad interp is good):
p = horzcat(s{:}); x = horzcat(p.x); nx = horzcat(p.nx);  % all nodes coords
dens = sin(1 + x(1,:) + 1.7*x(2,:) - 0.4*x(3,:))';   % any smooth surf func, col vec
k = 57;       % which pan to test on-surf u vals
nrs = [4 6 8 10];
for i=1:numel(nrs), o.nr = nrs(i); o.nt = 2*o.nr;
  s = add_panels_auxquad(s,o);                   % overwrites
  fprintf('setup interp... '); Linfo = setup_auxinterp(s{1}.t,o);  % slow
  u = LapDeval_panels(s{k},s,dens,Linfo);  % potential should be -1/2 at all pts
  fprintf('auxquad conv:   nr=%d\t   u(1) = %.16g\n',o.nr, u(1)) % should conv
end
% gets about 1e-10 at nr=8
