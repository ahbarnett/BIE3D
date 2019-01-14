function [f,fn,ft] = data_planewave(incdir,T,Tt,t,x,nx)
% DATA_PLANEWAVE  eval data f, fn, ft for plane wave solution to 3D wave eqn.
%
% [f,fn,ft] = data_planewave(incdir,T,Tt,t,x,nx)
%
% Inputs:
%  t - 1*n array of target times, or a single time
%  x, nx - 3*n arrays of target points and normal derivs (nx is optional)
%  incdir - 3-cmpt incident direction vector, need not be normalized.
%  T, Tt - source function T(s) of time, and (optionally) its deriv T'(s).
%
% Output:
%  f, fn, ft - n*1 col vecs of evaluated value, normal deriv at const t, and
%             t-deriv at const location.
%
% With no arguments, self-test done

% Barnett 1/13/19
if nargin==0, test_data_planewave; return; end

incdir=incdir(:)/norm(incdir);           % ensures wave travels at unit speed!
d = (incdir'*x)';
t = t(:);   % ensure col vec
f = T(t-d);
if nargout>1
  ft = Tt(t-d);
  fn = -(incdir'*nx)' .* ft; 
end
  
%%%%%%%%%
function test_data_planewave
T = @(t) t.^2; Tt = @(t) 2*t;                          % func of time, deriv
incdir = [2;1;2]/3;                                    % unit direction
x = [-1;2;0.5]; nx = [1;-2;2]/3; t = 1.2;              % arb targ, unit normal
[f,fn,ft] = data_planewave(incdir,T,Tt,t,x,nx);
'should be zero:'
f - T(t-dot(x,incdir))
ft - Tt(t-dot(x,incdir))
fn + dot(nx,incdir)*Tt(t-dot(x,incdir))
