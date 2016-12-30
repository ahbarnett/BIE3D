function a = relatedpanel(t,s)
% RELATEDPANEL  Identify if self, neighbor (give number), or far panel.
%
% a = relatedpanel(t,s)
% Inputs:
%  t : target panel struct (needs t.x and t.neix1, otherwise always ret 0)
%  s : source panel struct, or cell array of structs (each needs s.x)
% Output:
%  a : col vec with jth entry 0 if src pan j is far from t, or entry 1...9
%      describing which neighbor (2...9) or self (1) the jth src pan is of t.
%
% Relies on neix1, made in create_panels.m

% Barnett 12/29/12, improving on samepanel since no global pan indexing needed
n = numel(s);
a = zeros(n,1);
if ~isfield(t,'neix1'), return; end    % can't say anything about t's nei's
x = [t.x(:,1), t.neix1];    % 3-by-9 test list of coords for the one targ pan
if n==1 & ~iscell(s), s = {s}; end   % ensure cell array
for i=1:n
  j = find(dists(s{i}.x(:,1), x)<1e-14);
  if isempty(j), j=0; end
  a(i) = j;
end
