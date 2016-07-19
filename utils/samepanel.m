function a = samepanel(t,s)
% SAMEPANEL  Return true if two panels have numerically the same nodes
%
% a = samepanel(t,s)
% s may also be cell array of panels in which case a logical col vec

% Barnett 7/15/16
n = numel(s);
if n==1
  a = max(abs(s.x(1)-t.x(1)))<1e-14;   % cheapo test first pt only
else
  a = false(n,1);
  for i=1:n
    a(i) = max(abs(s{i}.x(1)-t.x(1)))<1e-14;   % cheapo test first pt only
  end
end
%a = (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14;
