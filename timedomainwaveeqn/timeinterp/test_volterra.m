% solve simple single-variable Volterra IE via simplified pred-corr variants.
% Barnett 12/23/16 based on Hagstrom description.
%
% We solve    u(t) + int_{t-1}^t u(s) ds = f(t),   for all t
% ie 2nd-kind Volterra w/ unit top-hat kernel from 0 to 1 unit into the past.
% We assume solution correct for t<0 then evolve for t>=0

clear
om = 2.0; soln = @(t) cos(om*t);    % true solution, freq om

% numerical params:
h = 0.1;        % timestep
m = 6;          % order of timestep scheme
p = 100;         % num Gauss-Legendre nodes to approx integral in the IE
% shouldn't dep so much on p! *** debug

fprintf('points per wavelength = %.3g\n',2*pi/om/h)
[y,w]=gauss(p); w = w/2; y = (y-1)/2;  % quadr scheme for smooth funcs on [-1,0]
[~,jmin,umat,upmat] = interpmat(y,h,m);  % jmin is most ancient grid offset
q = w*umat;     % add coeffs for each quadr node, scaled by its weight

nend = 100;   % how far to evolve starting at n=0
nn = jmin:nend;  tt = h*nn; % time indices, grid values
f = cos(om*tt) + (sin(om*tt)-sin(om*(tt-1)))/om;  % RHS driving eval on grid
utrue = soln(tt);
u = 0*nn;
known = (nn<0); %  indices for assumed known values, only t<0
u(known) = soln(tt(known));

% evolve some schemes... (note we solve for u at n=0,1,..)
disp('fully implicit..')    % no extrap needed
for n=0:nend, j = n+1-jmin;  % j is index in the u grid, n is the timestep
  u(j) = ( -sum(q(1:end-1).*u(j+jmin:j-1)) + f(j) )/(1+q(end)); % q(end)=q_0
end
fprintf('\tworst u err = %.3g\n',max(abs(u-utrue)))

%figure; subplot(2,1,1); plot(tt, [u;utrue;f], '.-'); legend('u','utrue','f'); subplot(2,1,2); semilogy(tt,abs(u-utrue),'.-'); legend('err');


