function A = intxperiinterpmat(N,Y,n,y,sphere)
% INTXPERIINTERPMAT Matrix that resamples 2d func on interval cross periodic.
%
% L = intxperiinterpmat(N,Y,n,y) returns the sum(N)-by-sum(n) dense matrix that
%  when acting on a smooth vector of values on a spectral grid on (x,y) in
%  [0,2pi)x[-1,1], returns the vector of interpolants on a new grid.
%  The input grid is described by y (the set of y values in [-1,1], on which a
%  "ring" of periodic points in x lies), and n (the corresponding numbers of
%  grid points per ring). The output grid is similarly given by Y and N.
%  The matrix has the same action as a call to intperiinterpmat.
%  See description of that latter function.
%
% L = intxperiinterpmat(N,Y,n,y,1) uses the variant where the y in [-1,1]
%  coordinate is assumed to be z=cos(theta) for a sphere parameterization,
%  so applies different rules to even and odd x-Fourier modes, as in
%  intxperiinterp.
%
% See also: intxperiinterp
%
% Note: dominant cost is the GEMMs below, ie blk B = LxI * F;  A = F * B;

% Alex Barnett 9/5/19. sphere sqrt(1-y^2) fix for odd modes 12/11/19.
if nargin==0, test_intxperiinterpmat; return; end
if nargin<5, sphere=0; end
if sum(mod(n,2)~=0)>0, error('all n must be even!'); end
if sum(mod(N,2)~=0)>0, error('all N must be even!'); end

nt = sum(n); Nt = sum(N);     % size of input and output spaces
M = max(N);                   % biggest output 1d-grid (ring) size
ny = numel(y); NY = numel(Y);
w = baryweights(y);  
L = baryprojs(y,w,Y);         % get 1d interp matrix Y <- y
LxI = kron(L,eye(M));         % (M*NY)-by-(M*ny), with 1:M fast on both axes
% Fill blk cols B by hit blk cols of LxI with small DFT matrices from right...
B = nan(M*NY,nt);
off = 0;
for i=1:ny                    % loop over input rings
  F = fft(eye(n(i))) / n(i);  % note quadr weight 1/n(i) for Euler-F
  if M>n(i)           % upsample, stacking to give M-by-n(i)
    F = [F(1:n(i)/2,:); F(n(i)/2+1,:)/2; zeros(M-n(i)-1,n(i)); F(n(i)/2+1,:)/2; F(n(i)/2+2:end,:)];
  else                % downsample or stay same (average makes symm, real->real)
    F = [F(1:M/2,:); (F(M/2+1,:)+F(end-M/2+1,:))/2; F(end-M/2+2:end,:)];
  end
  if sphere     % precorrect odd-m modes to make the func to interp in y smooth
    F = kron(ones(M/2,n(i)),[1; 1./sqrt(1-y(i).^2)]) .* F;
  end
  B(:,off+(1:n(i))) = LxI(:,M*(i-1)+(1:M)) * F;  % dominant cost, B takes vals -> F co
  off = off + n(i);
end
% Fill blk rows of A by hit blk rows of B with small DFT mats from left...
A = nan(Nt,nt);
off = 0;
for j=1:NY                    % loop over output rings, converting F coef to vals
  F = ifft(eye(N(j))) * N(j); % note weight N(i) for ifft
  if N(j)<M                   % downsample (average makes symm, real->real)
    F = [F(:,1:N(j)/2), F(:,N(j)/2+1)/2, zeros(N(j),M-N(j)-1), F(:,N(j)/2+1)/2, F(:,N(j)/2+2:end)];
  end
  if sphere     % postcorrect odd-m modes to recover their Y-func
    F = F .* kron(ones(N(j),M/2),[1, sqrt(1-Y(j).^2)]);
  end
  A(off+(1:N(j)),:) = F * B(M*(j-1)+(1:M),:);  % also dominant cost
  off = off + N(j);
end


%%%%%%
function test_intxperiinterpmat
ny = 30;                         % initial 1d grid
y = gauss(ny);                   % good interpolants for [-1,1]
n = ny*(2+y); n = ceil(n/2)*2;   % varying # in each row, make even
NY = 50;                         % final 1d grid
Y = gauss(NY);
N = NY*(2-Y); N = ceil(N/2)*2;   % a new variation in # in each row
for sphere=0:1                   % test both cases
  fprintf('case sphere = %d:\n',sphere)
  if sphere==0
    fun = @(x,y) exp(sin(x).*y); % 2pi-periodic in x, non-periodic in y, real
  else                           % trickier: test generic smooth func on S^2
    fR3 = @(a,b,c) 1./sqrt((2-a).^2+(1.5-b).^2+(-3-c).^2); % smooth R3 func 1/r
    fun = @(x,y) fR3(cos(x).*sqrt(1-y.^2), sin(x).*sqrt(1-y.^2), y); % restrict
  end
  f = intxperieval(fun,n,y);     % input data vector
  F = intxperieval(fun,N,Y);     % exact ans on final nodes
  tic %,profile clear; profile on
  A = intxperiinterpmat(N,Y,n,y,sphere);  % do it
  toc %,profile off; profile viewer
  g = A*f(:);                      % use it
  fprintf('A size %d-by-%d, max abs Im g=%d\n',size(A,1),size(A,2),max(abs(imag(g))))
  g = real(g);
  fprintf('interp rms err: %.3g\n', norm(g-F(:))/sqrt(numel(g)))
  %figure; plot(f); figure; plot(F); hold on; plot(g);  % debug

  tic; disp('compared vs crude feed-in-unit-vecs slower method...')
  A2 = nan(size(A));
  for j=1:sum(n)
    v = zeros(sum(n),1); v(j) = 1;
    A2(:,j) = intxperiinterp(v,N,Y,n,y,sphere);  % consistent against this code?
  end
  toc
  %figure; imagesc(imag(A2-A)); colorbar  % debug
  disp('errors vs the slower method:')
  max(abs(real(A2(:)-A(:)))), max(abs(imag(A2(:)-A(:)))), max(abs(A2*f(:)-g))
end
