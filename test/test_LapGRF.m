% test Laplace GRF in a torus for interior & singular on-surface cases
% Barnett 7/18/16

clear
% interior Laplace soln: cols in R3 mapped to values (row vec), and grad f...
lam = 1.0; f = @(x) exp(lam*x(1,:)).*cos(lam*x(2,:));
gradf = @(x) [lam*f(x); -lam*exp(lam*x(1,:)).*sin(lam*x(2,:)); 0*f(x)];
fprintf('did we get analytic grad f right? yes: %.3g\n',checkgrad(f,gradf))

% surface...
so.a=1; so.b = 0.5; [s N] = create_panels('torus',so,[]); % default # pans
[x nx] = getallnodes(s);
sigma = sum(nx.*gradf(x),1)'; tau = -f(x)';     % densities for GRF (col vecs)
showsurffunc(s,sigma); title('GRF: sigma');
showsurffunc(s,tau); title('GRF: tau');

% interior GRF test...
t.x = [0.9;-0.2;0.1]; t.N = 1;   % far targ pt (a fake panel), must be inside
u = LapSeval_panels(t,s,sigma) + LapDeval_panels(t,s,tau);  % S.sigma + D.tau
fprintf('N=%d, interior GRF test at 1 pt: u err = %.3g\n', N, u - f(t.x))
% about 1e-9

% on-surface GRF test using targets on one panel...
k = 57; t = s{k};            % the panel to use
o.nr = 8; o.nt = 2*o.nr;     % aux quad orders
s = add_panels_auxquad(s,o);
fprintf('setup interp... '); Linfo = setup_auxinterp(s{1}.t,o);
u = LapSeval_panels(t,s,sigma,Linfo) + LapDeval_panels(t,s,tau,Linfo);
fprintf('max on-surf GRF err on panel k = %.3g\n',max(abs(u - f(t.x)'/2)))
% about 1e-9
