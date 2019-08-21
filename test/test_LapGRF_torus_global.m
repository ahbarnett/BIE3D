% off-surface Green's rep formula test for any torus quadrature
% Barnett 8/16/19
clear

% interior Laplace soln: cols in R3 mapped to values (row vec), and grad f...
lam = 1.1; f = @(x) exp(lam*x(1,:)).*cos(lam*x(2,:));
gradf = @(x) [lam*f(x); -lam*exp(lam*x(1,:)).*sin(lam*x(2,:)); 0*f(x)];
fprintf('did we get analytic grad f right? yes: %.3g\n',checkgrad(f,gradf))

a=1.0; b=0.5;   % torus shape
% b = cruller(b,0.1,5,3); % optional, tougher test
xin = [0.9; -0.2; 0.1]; xout = [1.9; 0.7; 1.0];   % both "far" from surf
t.x = [xin,xout]; t.nx = randn(3,2);            % random target direc derivs
uex = f(t.x(:,1)); udex = t.nx(:,1)'*gradf(t.x(:,1));   % exact (val,dderiv)
for Na = 20:20:80, Nb = ceil(0.5*Na);      % tie minor discr to major
  s = setup_torus_doubleptr(a,b,[Na,Nb]);
  sigma = sum(s.nx.*gradf(s.x),1)'; tau = -f(s.x)';  % dens, int GRF (col vecs)
  [S Sd] = Lap3dSLPmat(t,s);
  [D Dd] = Lap3dDLPmat(t,s);
  u = S*sigma + D*tau;
  ud = Sd*sigma + Dd*tau;
  fprintf('N=[%d,%d] errs (val,grad):\tint %.3g,%.3g\text %.3g,%.3g\n',Na,Nb,u(1)-uex,ud(1)-udex,u(2),ud(2))   % int GRF gives zero outside
end
