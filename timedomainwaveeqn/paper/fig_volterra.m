% figure script for Volterra IE convergence for time-interp pred-corr scheme
% Barnett 11/14/18, just a m-loop (note m=2q) around test_volterra.m
%
% We solve    u(t) + int_{t-1}^t u(s) ds = f(t),   for all t
% ie 2nd-kind Volterra w/ unit top-hat kernel from 0 to 1 unit into the past.
% f(t) is built to match to the Gauss quadr accuracy the desired solution u(t).
% We assume solution correct for t<0 then evolve for t>=0. This allows u to
% be smooth, a requirement of the D-spline interp scheme.
% Note we don't shuffle the solution; we write once into a long array.
%
% We learn order is empirically 2q+2 = m+2, and that only the expected # p
% Gauss nodes are needed, although when p is small (10, sufficient for om=5),
% there is random error prefactor behavior but consistent with convergence ord.

clear
%om = 3.0; soln = @(t) cos(om*t);    % soln u(t) fun, const freq om (eg 3, 10)
om = 4.0; soln = @(t) exp(-(t-6).^2).*cos(om*t);  % wavepacket (starts 1e-16)
%solnp = @(t) -2*(t-6).*soln(t) - exp(-(t-6).^2).*sin(om*t);  % unused

Tend = 10.0;                        % when to stop, t

p = 16;             % fixed # Gauss-Legendre nodes to approx integral in the IE
[y,w]=gauss(p); w = w/2; y = (y-1)/2;  % quadr scheme for smooth funcs on [-1,0]

ms= 2:2:8;           % m=2q, orders of timestep scheme (really ord = m+2)
hs = 10.^-(0.7:.1:2);   % set of timestep sizes for convergence study

for k=1:numel(ms), m=ms(k)   % =========================== m loop
  for i=1:numel(hs), h = hs(i);      % -------------------- h convergence loop
    fprintf('h=%g  (points per wavelength = %.3g) :\n',h,2*pi/om/h)
    [~,jmin,umat,upmat] = interpmat(y,h,m);  % jmin is most ancient grid offset
    q = w*umat;     % add coeffs for each quadr node, scaled by its weight, row
    fprintf('jmin=%d, size(q) = %d\n',jmin,size(q))
    % note: these q are same as Tom's buildcofs, excluding the ident. u(t) term
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
    errs(k,i) = max(abs(u-utrue));
    fprintf('fully impl: \tu(T)=%.15g\tworst u err = %.3g\n',u(end),errs(k,i))
  end
end

% soln plot...
%h=0.1; nend = round(Tend/h); tt = h*nn; f = rhs(soln,tt); utrue = soln(tt);
figure; plot(tt, [utrue;f], '.-','markersize',3);
legend({'$u$','$f$'},'Interpreter','latex');
axis([0 Tend -1 1]); text(0.2,0.9,'(a)');
xlabel('$t$','interpreter','latex'); set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 volterra_soln.eps
%figure; semilogy(tt,abs(u-utrue),'.-'); legend('err');

% err conv plot...
figure; loglog(hs, errs, '+-'); axis tight; h0=0.4;  % h0 sets where laws hit 1
set(gca,'ColorOrderIndex',1); ht = hs(hs<0.1); % trunc list of hs for powers
hold on; loglog(ht, (ht/h0).^(ms(:)+2),'--');  % power laws
c1=num2cellstr(ms,[],'$2q=$ '); c2=num2cellstr(ms+2,2,'order ');
h=legend({c1{:},c2{:}});  % concatenate legend cell string list
set(h,'location','southeast','Interpreter','latex');
ylabel('max error $\|\tilde u -u \|_\infty$','interpreter','latex');
xlabel('$\Delta t$','interpreter','latex');
text(0.013,1e-3,'(b)');
set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 volterra_conv.eps

% plot of the q vector, showing D-splines totally unresolved by p-node rule...
figure; plot(q); title('q vector at the min dt (note: unresolved but good)');


%%%%%%
function f = rhs(ufunc, t)
% evaluate f(t), the RHS func, at an array of t values, to num quadr accuracy.
% See above defn of IE.
% Barnett 12/27/16. Option for lower limit max(0,t-1), 11/14/18.

%if numel(ufunc)>1,    % cell 2-array interpret as value and deriv.
%  ufuncp = ufunc{2};
%  ufunc = ufunc{1};
%else
%  ufuncp = @(t) 0*t;     % dummy deriv as zero func
%end

p=50;   % probably should differ from p in timestepping scheme, else crime.
[y,w]=gauss(p); w = w/2; y = (y-1)/2;  % quadr scheme for smooth funcs on [-1,0]
f = 0*t;
for i=1:numel(f)
  sc = 1; %max(0,min(1,t(i)));  % scale factor for interval (clips lower lim to 0)
  f(i) = ufunc(t(i)) + ...   % ufuncp(t(i)) +...      % for u + u' term, unused
         sc*w * ufunc(t(i)+sc*y);  % 2nd kind Volterra, do int
end
end
