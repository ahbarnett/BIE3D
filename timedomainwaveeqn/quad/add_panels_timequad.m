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

% Barnett 12/18/16, broke out disteval 12/29/16

for k=1:numel(pan)
  p = pan{k};
  if ~isfield(p,'naux'), error(sprintf('pan %d has no aux nodes!',k)); end
  n = p.naux;
  pan{k}.auxdelays = 0*p.auxwei;
  for j=1:p.N        % loop over targs in this panel
    pan{k}.auxdelays(:,j) = dists(p.auxnodes(:,:,j),p.x(:,j)); % col vec
  end
end
