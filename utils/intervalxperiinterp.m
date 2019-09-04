function g = intervalxperiinterp(f,N)
% INTERVALXPERIINTERP  resample 2d func on interval cross periodic to new grid.
%
% g = intervalxperiinterp(f,N,X,n,x) resamples f given on [-1,1]x[0,2pi)
%  at the 2d nodes of the following form:
%   . . . . . . . . . 
%   .  .  .  .  .  .  
%   ..................
%   |0                |2pi       <- this indicates horiz axis values
%  where the f values are ordered fast along each row, where the ith row
%  is a 0-indexed periodic grid for [0,2pi) with n(i) points, and the vertical
%  coordinate of all nodes in the ith row is X(i).
%  The output is a list of values g, corresponding to a grid of similar form,
%  given by lists N and X.
%
% Note on phasing: the output and input grid first entries align (ie, as if they
%  are both 0-indexed; note this matches setupquad), in the horizontal dimension.
%
% See also: BIE3D/utils/peri2dspecinterp

% Barnett 9/4/19
if nargin==0, test_intervalxperiinterp; return; end

n = size(f);                  % with luck a 2D array
if length(n)~=2, error('f must be a 2D array!'); end
if N==n, g = f; return; end   % trivial case!
if mod(N(1),2)~=0 || mod(N(2),2)~=0 || mod(n(1),2)~=0 || mod(n(2),2)~=0, warning('The two dims of both N and n must be even; sorry'); end
F = fft2(f);
% dim 1: pad or trunc? (note: unlike in 1D case, can also leave dim as is)
if N(1)>n(1)        % upsample
  F = [F(1:n(1)/2,:); F(n(1)/2+1,:)/2; zeros(N(1)-n(1)-1,n(2)); F(n(1)/2+1,:)/2; F(n(1)/2+2:end,:)];
else                % downsample or stay same
  F = [F(1:N(1)/2,:); F(end-N(1)/2+1:end,:)];
end
% dim 2: pad or trunc? (note: unlike in 1D case, can also leave dim as is)
if N(2)>n(2)        % upsample
  F = [F(:,1:n(2)/2) F(:,n(2)/2+1)/2 zeros(N(1),N(2)-n(2)-1) F(:,n(2)/2+1)/2 F(:,n(2)/2+2:end)];
else                % downsample or stay same
  F = [F(:,1:N(2)/2) F(:,end-N(2)/2+1:end)];
end
g = ifft2(F)*prod(N./n);   % go back, factors from the ifft


%%%%%%
function test_intervalxperiinterp
n = [60 100];          % starting sample rect grid

x = 2*pi*(0:n(1)-1)/n(1); y = 2*pi*(0:n(2)-1)/n(2);
[xx yy] = ndgrid(x,y);
f = @(x,y) exp(sin(x+2*y));   % doubly-periodic func
% some combos of upsample and downsample each dim...
Ns = {[98 100],[32 70],[32 144],[102 144],[60 130],[60 72]};
for i=1:numel(Ns), N=Ns{i};
  g = peri2dspecinterp(f(xx,yy),N);
  [xe ye] = ndgrid(2*pi*(0:N(1)-1)/N(1),2*pi*(0:N(2)-1)/N(2));
  ge = f(xe,ye);
  disp(norm(g(:) - ge(:)))
end
