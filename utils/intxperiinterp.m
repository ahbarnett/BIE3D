function g = intxperiinterp(f,N,Y,n,y,sphere)
% INTXPERIINTERP  resample 2d func on interval cross periodic to new grid.
%
% g = intxperiinterp(f,N,Y,n,y) spectrally resamples smooth f given on (x,y) in
%  [0,2pi)x[-1,1] at the 2d nodes of the following form:
% y=+1
%       . . . . . . . . .              <- n(3) pts on line y=y(3)
%       .  .  .  .  .  .               <- n(2) pts on line y=y(2)
%       ..................             <- n(1) pts on line y=y(1)
% y=-1
%       |x=0              |x=2pi
%  where the f values are ordered fast along each row, where the ith row
%  is a 0-indexed periodic grid for [0,2pi) with n(i) points, and the vertical
%  coordinate of all nodes in the ith row is -1 <= y(i) <= 1.
%  The output is a row vec of values g interpolated to a grid of similar form
%  given by lists N and Y. It will only be efficient if there is not too much
%  variation in the lists N or n. All entries of n and N are assumed even.
%
%  Note on phasing: output & input peri grid first entries align (ie, as if they
%  are both 0-indexed; note this matches setupquad), in x dimension.
%
% g = intxperiinterp(f,N,Y,n,y,1) instead assumes that the y ordinate is the
%  z=cos(theta) variable in a sphere parameterization, while x=phi is the
%  azimuthal variable, and that f(theta,phi) is smooth on the sphere. This
%  invokes a different interpolation rule: the even x-Fourier modes are
%  treated as smooth in y, whereas the odd x-Fourier modes are smooth when
%  divided by sqrt(1-y^2).
%  
%  Naming: "intxperi" = interval cross periodic = [-1,1] x [0,2pi)  = Y x X
%
% See also: BIE3D/utils/peri2dspecinterp, intxperieval, intxperiinterpmat

% Alex Barnett 9/4/19, symmetrized for reality 9/5/19.
% sphere sqrt(1-y^2) factor fix for odd modes 12/11/19.
if nargin==0, test_intxperiinterp; return; end
if nargin<6, sphere=0; end
if sum(mod(n,2)~=0)>0, error('all n must be even!'); end
if sum(mod(N,2)~=0)>0, error('all N must be even!'); end
f = f(:).';           % make row vec

% STEP 1: use 1d FFTs to fill coeffs array of sufficient x-freq-extent (height)
M = max(N);           % biggest output 1d-grid (ring) size, assume even
ny = numel(y);
fhat = zeros(M,ny);   % DFT coeffs, each col a ring, so contiguous 1d FFT access
off = 0;
for i=1:ny
  fihat = fft(f(off+(1:n(i)))) / n(i);  % row, note quadr wei 1/n(i) for Euler-F
  if M>n(i)           % upsample
    fhat(:,i) = [fihat(1:n(i)/2), fihat(n(i)/2+1)/2, zeros(1,M-n(i)-1), fihat(n(i)/2+1)/2, fihat(n(i)/2+2:end)];
  else                % downsample or stay same (average makes symm, real->real)
    fhat(:,i) = [fihat(1:M/2), (fihat(M/2+1)+fihat(end-M/2+1))/2, fihat(end-M/2+2:end)];
  end
  off = off + n(i);
end
% STEP 2: interpolate the DFT coeffs array in y, for each k_x coeff
w = baryweights(y);  
L = baryprojs(y,w,Y);         % get 1d interp matrix for smooth funcs Y <- y
if sphere                     % precorrect the odd-m modes to make smooth in y
  fhat = fhat .* kron(ones(M/2,1),[ones(1,ny); 1./sqrt(1-y(:)'.^2)]);
end
ghat = fhat*L.';              % interp each row of fhat to give corresp row ghat
NY = numel(Y);
if sphere                     % postcorrect the odd-m modes on the new Y-grid
  ghat = ghat .* kron(ones(M/2,1),[ones(1,NY); sqrt(1-Y(:)'.^2)]);
end
%figure; imagesc(log10(abs(fhat))); figure; imagesc(log10(abs(ghat))); % debug
clear fhat
% STEP 3: use 1d iFFTs to eval on each new ring
g = nan(1,sum(N));     % output list
off = 0;
for i=1:NY
  k = N(i)/2;          % half-size for output grid (size 2k ifft)
  % col vec of F coeffs len N(i), symmetrizes so preserves reality...
  co = [ghat(1:k,i); (ghat(k+1,i)+ghat(end-k+1,i))/2; ghat(end-k+2:end,i)];
  g(off+(1:N(i))) = ifft(co) * N(i);     % weight N(i)
  off = off + N(i);
end

%%%%%%
function test_intxperiinterp
ny = 20;                         % initial 1d grid
y = gauss(ny);                   % good interpolants for [-1,1]
n = ny*(2+y); n = ceil(n/2)*2;   % varying # in each row, make even
NY = 30;                         % final 1d grid
Y = gauss(NY);
N = NY*(2-Y); N = ceil(N/2)*2;   % a new variation in # in each row
for sphere=0:1                   % test both cases
  fprintf('case sphere = %d:\n',sphere)
  if sphere==0
    fun = @(x,y) exp(sin(x).*y); % 2pi-periodic in x, non-periodic in y, real
  else                           % trickier: test generic smooth func on S^2
    fR3 = @(a,b,c) a;            % simple smooth real func of R3
    fR3 = @(a,b,c) 1./sqrt((2-a).^2+(1.5-b).^2+(-3-c).^2); % smooth R3 func 1/r
    fun = @(x,y) fR3(cos(x).*sqrt(1-y.^2), sin(x).*sqrt(1-y.^2), y); % restrict
  end
  f = intxperieval(fun,n,y);     % input data vector
  F = intxperieval(fun,N,Y);     % exact ans on final nodes
  g = intxperiinterp(f,N,Y,n,y,sphere);   % do it
  fprintf('check reality: max abs Im g=%d\n',max(abs(imag(g))))
  g = real(g);
  fprintf('interp rms err: %.3g\n', norm(g-F)/sqrt(numel(g)))
end
figure; plot(f); figure; plot(F); hold on; plot(g);  % debug, crude 1d plots
