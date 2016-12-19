function pan = add_panels_timequad(pan, o)
% pan = add_panels_timequad(pan, o) adds time-domain info for each target in
%  each panel. Requires that add_panels_auxquad already run for all panels.
%
% Inputs:
%   pan = cell array of panels.
%   o = optional opts struct ... tdb
%
% Adds fields to each panel pan{k} in the output struct:
%   auxdelays (naux by targ) - time delays of each aux node for each target

% Barnett 12/18/16

for k=1:numel(pan)
  p = pan{k};
  if ~isfield(p,'naux'), error(sprintf('pan %d has no aux nodes!',k)); end
  n = p.naux;
  pan{k}.auxdelays = 0*p.auxwei;
  for j=1:p.N        % loop over targs in this panel
    pan{k}.auxdelays(:,j) = disteval(p.auxnodes(:,:,j),p.x(:,j)); % col vec
  end
end

function r = disteval(x,y)
% Dumb dense distance evaluator from sources to targets in R3.
% Inputs: x targ (3-by-M), y src (3-by-N)
% NB: fills M*N matrix, so assumes small cases
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
r = sqrt(d1.^2+d2.^2+d3.^2);          % matrix, size M*N