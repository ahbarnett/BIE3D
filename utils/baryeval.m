function u = baryeval(x, y, t)
%  u = baryeval(x, y, t) evaluates at the list of values t, the barycentric
%   interpolation formula using nodes x and with data y=f(x). Algorithm is O(NM)
%
% Based upon Berrut-Trefethen SIREV 2004, formula (4.2), 2nd true barycentric
%
% See also: BARYWEIGHTS, BARYPROJS

% Copyright (C) 2015 Alex Barnett

w = baryweights(x);
L = baryprojs(x, w, t); % matrix
u = L * y(:);                 % gives a col vec
u = reshape(u, size(t));
