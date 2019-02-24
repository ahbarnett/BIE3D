% test Laplace GRF in a torus for interior & singular on-surface cases
% Barnett 7/18/16. Included cruller, 2/24/19.

clear
% interior Laplace soln: cols in R3 mapped to values (row vec), and grad f...
lam = 1.0; f = @(x) exp(lam*x(1,:)).*cos(lam*x(2,:));
gradf = @(x) [lam*f(x); -lam*exp(lam*x(1,:)).*sin(lam*x(2,:)); 0*f(x)];
fprintf('did we get analytic grad f right? yes: %.3g\n',checkgrad(f,gradf))

% surface...
wobbly = 1;
so.a = 1; b = 0.5;   % major radius and poloidal radius
o.p = 8;   % panel order
if ~wobbly  % plain torus
  so.b =b;
else         % cruller
  wc = 0.1;  % surf modulation ampl
  wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
  wn = 5;   % # wobbles in toroidal, major
  g = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
  gt = @(t,p) -wc*wm*sin(wm*t+wn*p); gp = @(t,p) -wc*wn*sin(wm*t+wn*p);
  so.b = {g,gt,gp};     % pass in instead of b param
end
%so.np = 12; so.mp = 8;        % set # pans in each direction. default
so.np = 18; so.mp = 12;        % set # pans in each direction
[s N] = create_panels('torus',so,o);
[x nx] = getallnodes(s);
sigma = sum(nx.*gradf(x),1)'; tau = -f(x)';     % densities for GRF (col vecs)

% interior GRF test...
t.x = [0.9;-0.2;0.1]; t.N = 1;   % far targ pt (a fake panel), must be inside
u = LapSeval_panels(t,s,sigma) + LapDeval_panels(t,s,tau);  % S.sigma + D.tau
fprintf('N=%d, interior GRF test at 1 pt: u err = %.3g\n', N, u - f(t.x))
% about 1e-9

showsurffunc(s,sigma); title('GRF: sigma'); hold on; plot3(t.x(1),t.x(2),t.x(3),'k.');  % include interior targ pt
showsurffunc(s,tau); title('GRF: tau'); hold on; plot3(t.x(1),t.x(2),t.x(3),'k.');

% on-surface GRF test using targets on one panel...
k = 57; t = s{k};            % the panel to use
o.nr = 2*o.p; if wobbly, o.nr=2*o.p+4; end     % radial aux quad order
o.nt = 2*o.nr;     % angular aux quad order
s = add_panels_auxquad(s,o);
fprintf('setup interp... '); Linfo = setup_auxinterp(s{1}.t,o);
tic; u = LapSeval_panels(t,s,sigma,Linfo) + LapDeval_panels(t,s,tau,Linfo);
t0=toc;
fprintf('1-panel time %.3g s (%.3g targs/sec)\n',t0,o.p^2/t0)
fprintf('max on-surf GRF abs err on panel k = %.3g\n',max(abs(u - f(t.x)'/2)))
% about 1e-9 for torus

% whole-surface GRF...
tic;
u = LapSeval_panels(s,s,sigma,Linfo) + LapDeval_panels(s,s,tau,Linfo);
t0=toc;
fprintf('max on-surf GRF rel err = %.3g\n',max(abs(u - f(x)'/2))/max(abs(f(x))))
% torus p=8,np=12,mp=8: about 2sec, 4e-8
fprintf('time %.3g s (%.3g targs/sec)\n',t0,N/t0)

% cruller: roughly 1e3 targs/sec @ 1e-6 rel err

