function u = RetDeval_panels(tpan,span,dens,paninfo)
% RETDEVAL_PANELS  evaluate extra retarded wave eqn DLP term from surface in R3
% NOT YET incl on-surface
%
% u = RetDeval_panels(t,s,dens,paninfo)
%
% Called with no arguments does a self-test (see bottom of code for usage)
%
% dens is time deriv of retarded density tau_t(t-|x-y|,y), as func of y on Gamma
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

% Barnett 12/9/16 based on LapDeval_panels
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
    u(tinds) = RetDLPeval(t.x,y,ny,w,dens);
  else                          % one src panel is self (same as targ)
    ss = span{spanselfk};       % the self src pan
    far = true(Np,1); far(ss.nei) = false; far(spanselfk) = false;
    [y ny w] = getallnodes({span{far}});    % note makes back a cell array
    farinds = logical(kron(far, ones(n,1)));  % assumes same n (nodes) for all src pans
    u(tinds) = RetDLPeval(t.x,y,ny,w,dens(farinds));
    nearrd = dens(~farinds);
    auxdens = zeros(1,ss.naux*t.N);  % step 1: interp dens to all aux nodes:
    for i=1:9                   % loop over self+nei 3x3 local (near) pans
      if i==1, k=spanselfk; else, k=span{spanselfk}.nei(i-1); end % pan index
      densk = dens(span{k}.inds);   % get density data for this pan
      auxdens(paninfo.auxinds{i}) = paninfo.L{i} * densk;
    end
    auxdens = reshape(auxdens,[ss.naux, t.N]);  % make easier to index
    unear = nan(t.N,1);            % step 2: apply kernel at aux nodes
    for i=1:t.N             % loop over targs; each has own set of aux nodes
      unear(i) = RetDLPeval(t.x(:,i),squeeze(ss.auxnodes(:,:,i)),squeeze(ss.auxnormals(:,:,i)),ss.auxwei(:,i)',auxdens(:,i));
    end
    u(tinds) = u(tinds) + unear;
  end
  indoff = indoff + t.N;
end

function u = RetDLPeval(x,y,ny,w,dens)
% Dumb dense Retarded DLP extra term evaluator from sources density tau_t
% w/ given weights, to targets
% Inputs:
% x targ (3-by-M), y src, ny src normals (both 3-by-N), w src weights (1-by-N)
%
% NB: fills M*N matrix, so assumes small cases
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
rr = d1.^2+d2.^2+d3.^2;   % M*N
ddotn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
A =  ddotn ./ rr;  % kernel matrix w/o 1/4pi or src quad weights: cos(theta)/r
u = (1/4/pi) * (A * (w(:).*dens));
