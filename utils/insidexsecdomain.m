function i = insidexsecdomain(shape,so,x,z)
% INSIDEXSECDOMAIN  hacky routine to return boolean if (x,z) pts inside slice
%
% i = insidexsecdomain(shape,so,x,z)
%  where shape and so are the surface params as in create_panels, and x,z are
%  equal-shaped lists of coordinates, then returns boolean i, same shape as
%  x or z, which is true if i is inside the surface, false if outside.

% Hacky.  Barnett 1/13/19

if size(x)~=size(z), error('sizes of x and z must match!'); end
i = false(size(x));
if strcmp(shape,'torus')
  th = 0*x;     % thetas for all pts
  for s = [-1,1]        % handle two sides (signs) separately
    ii=(s*x>=0);        % indices of output pts on this side
    th(ii) = atan2(z(ii),s*x(ii)-so.a);
    th(isnan(th)) = 0;
    phi=angle(s);
    thunr = th(ii(:)');   % unrolled thetas for this size, row vec
    xx = torusparam(so.a,so.b,phi*ones(size(thunr)),thunr);   % get true pts
    i(ii) = dists(xx,[s*so.a;0;0]) > sqrt((x(ii(:))-s*so.a).^2+z(ii(:)).^2);
  end
else
  error('unknown shape!');
end
