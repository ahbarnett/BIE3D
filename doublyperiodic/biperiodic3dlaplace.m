function biperiodic3dlaplace
% Evaluate the potential due to N monopoles or dipoles with bi-periodic Laplace
% kernel in 3D, & check periodicity. No quadrature yet. No FMM yet.
% Barnett 13/3/1 - 13/3/5. collab Wilkening.
%testDLPderiv, return

N = 1e3; v = 2;  % problem size # source pts; verbosity (0,1,2,...)
fmm = 0; % use FMM to apply (1) or direct N^2 summation (0)
srctype = 'd'; % source type: 's' for monopoles (SLP), 'd' for dipoles (DLP)
s.x = rand(3,N)-kron([.5;.5;.5],ones(1,N)); % rand src pts in box
s.x(3,:) = s.x(3,:)/2; % keep max z away from zU and zD
%N=2; s.x = [.5 .4;.5 .4;0.25 0.25]; % ...or small test src pair
s.nx = rand(3,N)-kron([.5;.5;.5],ones(1,N)); % random normals (only used by DLP)
s.nx = s.nx .* repmat(1./sqrt(sum(s.nx.^2,1)),[3 1]); % unit length
coef = 2*(rand(N,1)-0.5); % source strength coeffs col vec
if srctype=='s', if mod(N,2), warning('N odd so can''t balance monopoles'); end
  coef(1+N/2:end) = -coef(1:N/2); end % enforce charge balance for monopoles  
  
d=1; e1 = [d;0;0]; e2 = [0;d;0]; % unit cell square for now, d-by-d
zU = 0.6; zD = -0.6; % box top and bottom face heights (larger makes k conv)
M = 22;   % # gauss nodes in z direction on faces L,R,etc
Mu = 20;  % gauss nodes per line segment in xy plane
Mp = 22;  % periodic nodes per edge on U,D faces (should be > 2K)
K = 10;   % number of Fourier terms either side of zero in each direc
nei = 1;  % # directly summed neighbors each side (0,1,2; 1 is usual 3x3)

[L R B T U D Mf] = buildbox(M,Mu,Mp,zU,zD,d); % 6 faces for collocation on
p.x = sphpts(24, 1.4*d); % MFS pts for representing periodizing part of field
if v,figure; showface({D U L R B T}); hold on; plot3(s.x(1,:),s.x(2,:),s.x(3,:),'k.');plot3(p.x(1,:),p.x(2,:),p.x(3,:),'r.'); axis tight; end

tic, [AL ALn] = SLPmatrix(L,p); [AR ARn] = SLPmatrix(R,p); % periodize: build Q
[AB ABn] = SLPmatrix(B,p); [AT ATn] = SLPmatrix(T,p);
Q = [AL-AR; ALn-ARn; AB-AT; ABn-ATn]; clear AL ALn AR ARn AB ABn AT ATn;
[AU AUn] = SLPmatririx(D,p);
V = [AU; AUn; AD; ADn]; clear AU AUn AD ADn; % same ordering as for W matrix
[kx ky] = meshgrid(-K:K,-K:K); krb = (2*pi/d)*(kx(:)+1i*ky(:)); % C plane k's
kap = sqrt(-abs(krb).^2); Nb = numel(krb); % radiative modes z-wavenumbers
W = zeros(size(V,1), 2*Nb); Mpp = Mp*Mp; % nodes per up face
for j=1:Nb        % fill each column at once, stack order: U Un D Dn
  kx=real(krb(j)); ky=imag(krb(j)); 
  W(1:Mpp,j) = exp(1i*([kx,ky,0]*U.x)'); % note z-zU=0 on entire face
  W(Mpp+1:2*Mpp,j) = W(1:Mpp,j) .* (1i*kap(j)); % up-normal derivs (radiate up)
end  % probably should recode to make all real numbers, 3x faster
W(2*Mpp+1:3*Mpp,1+Nb:end) = W(1:Mpp,1:Nb); % since D same as U
W(3*Mpp+1:end,1+Nb:end) = -W(1+Mpp:2*Mpp,1:Nb); % - since rad down but normal up
E = [Q, zeros(size(Q,1), size(W,2)); V, -W];  % stack periodizing matrix
fprintf('E size %dx%d, filled in %.3g sec...\n', size(E,1),size(E,2),toc)
if v>2,figure; imagesc(log10(abs(E))); colorbar; title('E matrix, log10'); end
tic; [EU ES EV] = svd(E,'econ'); fprintf('svd(E) in %.3g sec...\n',toc)
% figure; imagesc(abs(EV)); colorbar;
ESU = inv(ES)*EU'; clear ES EU; % now EV*(ESU* . ) does bkwds-stable pinv(E)

% Do periodizing of the N srcs:
tic, rhs=zeros(size(E,1),1); xo = s.x; % save original src locations
for i=-nei:nei, for j=-nei:nei % cancelling numerical for now (analytic better!)
    s.x = xo + repmat(i*e1 + j*e2, [1 N]);  % overwrite src w/ direct near copy
    if srctype=='s', [uL uLn] = SLPmatrix(L,s); [uR uRn] = SLPmatrix(R,s);
      [uB uBn] = SLPmatrix(B,s); [uT uTn] = SLPmatrix(T,s);
    else [uL uLn] = DLPmatrix(L,s); [uR uRn] = DLPmatrix(R,s);
      [uB uBn] = DLPmatrix(B,s); [uT uTn] = DLPmatrix(T,s); end
    discrep = [uL-uR; uLn-uRn; uB-uT; uBn-uTn] * coef; % same as in Q
     if srctype=='s', [uU uUn] = SLPmatrix(U,s); [uD uDn] = SLPmatrix(D,s);
       else [uU uUn] = DLPmatrix(U,s); [uD uDn] = DLPmatrix(D,s); end
    mismTB = [uU; uUn; uD; uDn] * coef;
    rhs = rhs - [discrep; mismTB]; % note minus sign since 
end, end, fprintf('box data eval in %.3g sec\n', toc)
cop = EV*(ESU*rhs);  % solve for periodizing degrees of freedom
fprintf('cop norm = %.3g, resid norm = %.3g\n', norm(cop), norm(E*cop-rhs))
if v>1,figure; semilogy(abs(cop)); title('coeffs: MFS then RBs'); end
coa = cop(1:size(Q,2)); corb = cop(size(Q,2)+1:end); % MFS and up-down RB coeffs
% now eval and check total field is periodic at points t nr corner...
x = [-d/2+0.05; -d/2-0.03; 0.1]; t.x = [x, x+e1, x+e2];
u = SLPmatrix(t,p) * coa; % periodizing part of field
for i=-nei:nei, for j=-nei:nei  % plus sum over direct near copies
  s.x = xo + repmat(i*e1 + j*e2, [1 N]);
  if srctype=='s', u = u + SLPmatrix(t,s) * coef;
  else u = u + DLPmatrix(t,s) * coef; end
  end, end
fprintf('potential periodic err, 2 direcs: %.3g, %.3g\n', u(2)-u(1), u(3)-u(1))

if v % Evaluate on a xy-plane slice at height z0:
nx = 100; z0 = -0.4; g = d*((0.5:nx-0.5)/nx-0.5);
[xx yy] = meshgrid(g,g); t.x = [xx(:)'; yy(:)'; z0+0*xx(:)']; % targets
tic, ug = SLPmatrix(t,p) * coa; % slice potentials, periodizing part
for i=-nei:nei, for j=-nei:nei  % add sum over direct near copies
    s.x = xo + repmat(i*e1 + j*e2, [1 N]);
    if srctype=='s', ug = ug + SLPmatrix(t,s) * coef;
      else ug = ug + DLPmatrix(t,s) * coef; end
  end, end, s.x = xo; % reset srcs
fprintf('eval on %dx%d plane done in %.3g sec.\n',nx,nx,toc)
up = zeros(3*nx,3*nx); % periodize z-slice to check it's per
Xp = up; Yp = up; % fill a larger 3x3 array
for i=-1:1, for j=-1:1, ii = (1:nx)+(i+1)*nx; jj = (1:nx)+(j+1)*nx; % inds
    Xp(jj,ii) = xx+i*e1(1)+j*e2(1); Yp(jj,ii) = yy+i*e1(2)+j*e2(2);
    up(jj,ii) = reshape(ug,size(xx)); end, end  % note jj=x direc, ii=y direc
figure; plot3(s.x(1,:),s.x(2,:),s.x(3,:),'k.'); hold on; showface({L R T B});
surf(Xp,Yp,z0+0*Xp,real(up)); %v = max(caxis); caxis(v*[-1 1]);
colormap(jet(256)); colorbar; shading interp; axis equal tight vis3d;
end

keyboard % allow user to examine variables; use dbquit to exit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A An] = SLPmatrix(t,s)
% matrix mapping charge strengths to potentials (and optional targ n-deriv)
% s.x are src locations; t.x t.nx are targets and normals
M = size(t.x,2); N = size(s.x,2); d = zeros(3,M,N);    % displacements
for i=1:3, d(i,:,:) = repmat(t.x(i,:)',[1 N]) - repmat(s.x(i,:),[M 1]); end
r = reshape(sqrt(sum(d.^2,1)), [M N]);   % distance matrix
if nargout==1
  A = 1./((4*pi)*r);
elseif nargout>1               % want n-deriv too
  ir = 1./r; A = ir/(4*pi);
  tndotd = 0*r; % fill matrix of targ normals dot displacements...
  for i=1:3, tndotd=tndotd+reshape(d(i,:,:),[M N]).*repmat(t.nx(i,:)',[1 N]);end
  An = -tndotd .* ir.^2 .*A; % -(1/4pi)(x.n)/r^3
end

function [A An] = DLPmatrix(t,s)
% matrix mapping dipole strengths to potentials (and optional n-deriv)
% s.x s.nx are src dipole locs and direcs; t.x t.nx are targets and normals
M = size(t.x,2); N = size(s.x,2); d = zeros(3,M,N);    % displacements
for i=1:3, d(i,:,:) = repmat(t.x(i,:)',[1 N]) - repmat(s.x(i,:),[M 1]); end
r = reshape(sqrt(sum(d.^2,1)), [M N]);   % distance matrix
ir = 1./r;
sndotd = 0*r; % fill matrix of src normals dot displacements...
for i=1:3, sndotd = sndotd+reshape(d(i,:,:),[M N]).*repmat(s.nx(i,:),[M 1]); end
if nargout==1
  A = (1/4/pi) * ir.^3 .* sndotd; % +(1/4pi)(x.n)/r^3
elseif nargout>1               % want n-deriv too
  irr = ir.^2;
  temp = (1/4/pi) * ir.*irr;
  A = temp .* sndotd; % +(1/4pi)(x.n)/r^3
  tndotd = 0*r;
  for i=1:3, tndotd=tndotd+reshape(d(i,:,:),[M N]).*repmat(t.nx(i,:)',[1 N]);end
  tndotsn = 0*r;
  for i=1:3, tndotsn=tndotsn+repmat(t.nx(i,:)',[1 N]).*repmat(s.nx(i,:),[M 1]);end
  An = temp .* (-3*irr.*tndotd.*sndotd + tndotsn); % dipole deriv
end

function [L R B T U D Mf] = buildbox(M,Mu,Mp,zU,zD,d);
% build six face structs, and Mf = # nodes per side face
x = d * (((1:Mp)-0.5)/Mp - 0.5);   % uniform in [-.5,.5], for peri trap rule
[xx yy] = meshgrid(x,x);
U.x = [xx(:)'; yy(:)'; zU+0*xx(:)']; U.nx = [0*xx(:)'; 0*xx(:)'; 1+0*xx(:)'];
D.x = [xx(:)'; yy(:)'; zD+0*xx(:)']; D.nx = U.nx; % both normals point up
[x w] = gauss(M); wz = w*(zU-zD)/2;  % vertical quadr weights
z = zD + (zU-zD)*(x+1)/2;    % node z-coords on side faces
zz = repmat(z',[Mu 1]); Mf = numel(zz); % nodes per side face
x = d*gauss(Mu)'/2; % for L and B walls
L.x = [-ones(1,Mf)*d/2; kron(ones(1,M), x); zz(:)'];
L.nx = [ones(1,Mf); zeros(1,Mf); zeros(1,Mf)];
R = L; R.x = R.x+repmat([d;0;0], [1 Mf]);
B.x = [kron(ones(1,M), x); -ones(1,Mf)*d/2; zz(:)'];
B.nx = [zeros(1,Mf); ones(1,Mf); zeros(1,Mf)];
T = B; T.x = T.x+repmat([0;d;0], [1 Mf]);

function [x,w] = gauss(N) % Legendre nodes and weights on [-1,1] from Trefethen
beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D); [x,i] = sort(x);
w = 2*V(1,i).^2;

function showface(Flist)
% SHOWFACE - plot face (or cell array of faces) nodes & surface normals, in 3d
l = 0.1;  % normal plotting length
c = 'bgrcmyb';   % color ordering
for i=1:numel(Flist)
  if iscell(Flist), F = Flist{i}; else F = Flist(i); end
  plot3(F.x(1,:),F.x(2,:),F.x(3,:), 'k.'); hold on;
  plot3([F.x(1,:); F.x(1,:)+l*F.nx(1,:)], [F.x(2,:); F.x(2,:)+l*F.nx(2,:)], [F.x(3,:); F.x(3,:)+l*F.nx(3,:)], [c(i) '-']);
end
axis equal vis3d;

function x = sphpts(n,R) % points on sphere rad R at origin, w/ 2n pts on circ
z = gauss(n);
mmin = 10; % min number of nodes on a latitude line (dep on desired accuracy)
x = [];  % pts on unit sphere
for i=1:n, r = sqrt(1-z(i)^2); % this xy rad
  m = sqrt((2*n*r)^2 + mmin^2); t = (1:m)/m*2*pi; % azimuthal angles
  x = [x, [r*cos(t); r*sin(t); 0*t+z(i)]];
end
x = R*x; % scale

function testevalspeed
% test speed of eval prozy field directly, similar speed to locexp?
M = 1e3;
N = 1e4;
y = rand(3,M);
s.x = rand(3,N);
co = rand(N,1)+1i*rand(N,1);
u = zeros(N,1);
tic, for m=1:M
  u = u + co(m)*(1./sqrt(sum((s.x-repmat(y(:,m),[1 N])).^2,1)))';
end
toc

function testDLPderiv   % check dipole same as monopole finite diff pair
x = [1;2;3]; n = [sqrt(1/3);sqrt(2/3);sqrt(2/3)];
e = 1e-5; s.x = [x-e*n/2,x+e*n/2]; s.nx = [n,n]; t.x = [0;1;0]; t.nx = [1;0;0];
[A An] = SLPmatrix(t,s); co = [-1;1]/e; % approx a dipole by monopoles
s.x = x; s.nx = n; [A1 An1] = DLPmatrix(t,s);
A*co - A1, An*co - An1
