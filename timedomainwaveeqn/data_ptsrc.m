function [f,fn,ft] = data_ptsrc(xs,T,Tt,t,x,nx)
% DATA_PTSRC  eval data f, fn, ft for fixed-loc t-dep pt src, 3D wave equation.
%
% [f,fn,ft] = data_ptsrc(xs,T,Tt,t,x,nx)
%
% Inputs:
%  t - 1*n array of target times, or a single time
%  x, nx - 3*n arrays of target points and normal derivs (nx is optional)
%  xs - 3*1 source location, fixed
%  T, Tt - source function T(s) of time, and (optionally) its deriv T'(s).
%
% Output:
%  f, fn, ft - n*1 col vecs of evaluated value, normal deriv at const t, and
%             t-deriv at const location.
%
% With no arguments, self-test done

% Barnett 12/15/16
if nargin==0, test_data_ptsrc; return; end

t= t(:)';   % ensure row vec
d = bsxfun(@minus,x,xs);   % displacements of targs rel to src
r = sqrt(sum(d.^2,1));
ir = 1./r;
Tret = T(t-r);
Ttret = Tt(t-r);
f = (1/4/pi)*(ir.*Tret)';
if nargout>1
  cth = sum(bsxfun(@times,nx,d),1) .* ir;
  ft = (1/4/pi)*(ir.*Ttret)';
  fn = (-1/4/pi)*(cth.*(ir.*ir.*Tret + ir.*Ttret))';   % two terms
end

%%%%%%%%%
function test_data_ptsrc
T = @(t) t.^2; Tt = @(t) 2*t;                          % func of time, deriv
xs = [1;2;3];                                          % src
x = [2;2;3]; nx = [1;0;0]; t = 1.2;                    % targ, dist 1
[f,fn,ft] = data_ptsrc(xs,T,Tt,t,x,nx);                % recall wavespeed = 1
tret = t-1;
'should be zero:'
f-T(tret)/(4*pi)
fn-(-Tt(tret)-T(tret))/(4*pi)
ft-Tt(tret)/(4*pi)
