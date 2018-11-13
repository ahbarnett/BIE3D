function rho = maxeigpowermeth(Afun,N,its,pad)
% MAXEIGPOWERMETH.  Estimate largest eigenvalue magnitude of linear operator
%
% rho = maxeigpowermeth(Afun,N,its,pad) estimates largest eigenvalue
%  magnitude of N*N matrix applied by mat-vec function handle Afun.
% its iterations are used, and the slope estimates using the last pad
% of them.

% Barnett 11/13/18

if nargin==0, test_maxeigpowermeth; return; end

x = randn(N,1);
for i=1:its
  x = Afun(x);
  if i==its-pad
    n1 = norm(x);
  end
end
n2 = norm(x);
rho = (n2/n1).^(1/pad);  % mean growth rate over [its-pad,its]

%%%%%%%%%
function test_maxeigpowermeth
n = 100;

A = randn(n);  % small test...
rho = maxeigpowermeth(@(x) A*x,n,100,10);
d = max(abs(eig(A)));
fprintf('randn mat: max eig mag: true %.6g, est %.6g \t(rel err: %.3g)\n',d,rho,abs(d-rho)/d)

n = 1000;  % meatier test...
A = diag(ones(n-1,1),-1); A(1,:) = 0.01*randn(1,n);  % weak companion mat (weaker seems harder)
A(1,1) = 0.9;
A = sparse(A);
rho = maxeigpowermeth(@(x) A*x,n,300,100);
lam = eig(full(A)); figure; plot(lam,'+'); axis equal; drawnow
d = max(abs(lam));  % back to full
fprintf('companion mat: max eig mag: true %.6g, est %.6g \t(rel err: %.3g)\n',d,rho,abs(d-rho)/d)
darp = max(abs(eigs(A,30,'LM',struct('tol',1e-6,'disp',0)))); % biggest 6 in mag
fprintf('\t\t\tcompare ARPACK eigs: %.6g \t(rel err: %.3g)\n',darp,abs(d-darp)/d)
