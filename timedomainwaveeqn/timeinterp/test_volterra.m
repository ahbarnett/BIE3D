function test_volterra(verb)
% test_volterra(verb)
% verb = 0: text only, 1: check p-conv of Gauss & final fig, 2: more figs
%
% solve simple single-variable Volterra IE via simplified pred-corr variants,
% and Hagstrom spline weights (interpmat call).
%
% We solve    u(t) + int_{t-1}^t u(s) ds = f(t),   for all t
% ie 2nd-kind Volterra w/ unit top-hat kernel from 0 to 1 unit into the past.
% f(t) is built to match to the Gauss quadr accuracy the desired solution u(t).
% We assume solution correct for t<0 then evolve for t>=0.
% Note we don't shuffle the solution; we write once into a long array.
%
% We learn order is empirically m+2, and that only the expected # p Gauss nodes
% are needed, although when p is small (10, sufficient for om=5), there is
% random error prefactor behavior but consistent with convergence order.
%
% Barnett 12/23/16. 12/27/16 general soln case, pred-corr, plots

if nargin==0, verb = 0; end
om = 5.0; soln = @(t) cos(om*t);    % true solution func, freq om (eg 3, 10)
Tend = 10.0;                        % what time to stop

if verb
  disp('convergence study of quadr scheme used for integral in IE...')
  ps = (5:5:30)'; Is=nan*ps; for i=1:numel(ps)
    [y,w]=gauss(ps(i)); w=w/2; y=(y-1)/2; % scheme for smooth funcs on [-1,0]
    Is(i) = w*soln(y);
  end, [ps(1:end-1), abs(Is(1:end-1)-Is(end))]  % results table
end

% numerical params:
hs = 1./(5:5:40);   % set of timesteps delta.t for convergence study
m = 6;          % order of timestep scheme
p = 16;         % num Gauss-Legendre nodes to approx integral in the IE

[y,w]=gauss(p); w = w/2; y = (y-1)/2;  % quadr scheme for smooth funcs on [-1,0]
errs = nan(2,numel(hs));               % save for plot

for i=1:numel(hs), h = hs(i);      % -------------------- h convergence loop
  fprintf('h=%g  (points per wavelength = %.3g) :\n',h,2*pi/om/h)
  [~,jmin,umat,upmat] = interpmat(y,h,m);  % jmin is most ancient grid offset
  q = w*umat;     % add coeffs for each quadr node, scaled by its weight, row
  % note: these q are same as Tom's buildcofs, excluding the identity u(t) term
  nend = round(Tend/h);                 % how far to evolve starting at n=0
  nn = jmin:nend; tt = h*nn;            % time indices, grid values
  f = rhs(soln,tt); utrue = soln(tt);   % RHS and true soln on grid
  fprintf('u_ex(T)=%.15g\n',utrue(end))
  known = (nn<0);  % indices for assumed known values, only t<0

  % evolve some schemes:
  u = 0*nn; u(known) = soln(tt(known));    % FULLY IMPL: set up soln array
  for n=0:nend, j = n+1-jmin;  % j is index in the u grid, n is timestep index
    u(j) = ( -sum(q(1:end-1).*u(j+jmin:j-1)) + f(j) )/(1+q(end)); % q(end)=q_0
  end
  errs(1,i) = max(abs(u-utrue));
  fprintf('fully impl: \tu(T)   =%.15g\tworst u err = %.3g\n',u(end),errs(1,i))
  
  u = 0*nn; u(known) = soln(tt(known));    % PRED-CORR: set up soln array
  we = extrap(m);                          % Tom's extrap weights, row
  for n=0:nend, j = n+1-jmin;
    uextrap = we*u(j-numel(we):j-1)';
    u(j) = f(j) - q*[u(j+jmin:j-1), uextrap]';         % predictor
    u(j) = u(j) - q(end)*(u(j)-uextrap);               % corrector
  end
  errs(2,i) = max(abs(u-utrue));
  fprintf('pred-corr: \tu(T)   =%.15g\tworst u err = %.3g\n',u(end),errs(2,i))
  if verb>1, figure; subplot(2,1,1); plot(tt, [u;utrue;f], '.-'); legend('u','utrue','f'); title(sprintf('h=1/%d',1/h)); subplot(2,1,2); semilogy(tt,abs(u-utrue),'.-'); legend('err'); end
end                                % -------------------

if verb, figure; loglog(1./hs, errs, '+-'); xlabel('1/h');ylabel('max u err');
  hold on; fac = 0.5; plot(1./hs, (fac*om*hs).^(m+2), 'r-'); axis tight;
  legend('fully impl','pred-corr','order m+2');
end

%%%%%%
function f = rhs(ufunc, t)
% evaluate f(t), the RHS func, at an array of t values, to num quadr accuracy
% Barnett 12/27/16
p=50;   % probably should differ from p in timestepping scheme, else crime.
[y,w]=gauss(p); w = w/2; y = (y-1)/2;  % quadr scheme for smooth funcs on [-1,0]
f = 0*t;
for i=1:numel(f)
  f(i) = ufunc(t(i)) + w * ufunc(t(i)+y);    % 2nd kind Volterra, do integral
end
