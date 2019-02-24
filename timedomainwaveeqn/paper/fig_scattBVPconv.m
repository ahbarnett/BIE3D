% final c,d subfigs of Fig.9, signal and convergence of scatt BVP example
% Barnett 1/15/19

clear
%load('../expts/wobblytorus/scattBVPconv_npupto12.mat'); % by:gen_scattBVPconv.m
%load('../expts/wobblytorus/scattBVPconv_dt0.9h_nrbetter.mat'); % by:gen_scattBVPconv.m
%load('../expts/wobblytorus/scattBVPconv_wider_dt1.0h.mat'); % by:gen_scattBVPconv.m
load('../expts/wobblytorus/scattBVPconv_wider_big_dt1.0h.mat'); % by:gen_scattBVPconv.m
%load('../expts/wobblytorus/scattBVPconv_dt0.8h.mat'); % by:gen_scattBVPconv.m
% must match...
dx=0.05; gx=-2:dx:2; gz = -1:dx:2; [xx zz] = meshgrid(gx,gz);
t.x = [xx(:)';0*xx(:)';zz(:)']; t.N=numel(xx);           % x-z plane of targs
jx=37; jz=41; j=jz+numel(gz)*(jx-1); x0=xx(j); z0=zz(j);  % spatial pt
fprintf('test pt x0=(%g,%g,%g)\n',x0,0,z0);

% ================ (c)
r = run{end};  % run to show signal
N = size(r.muall,1);
figure(1); clf; hold off; sig = r.utot(jz+numel(gz)*(jx-1),:);  % signal at x0
%plot(r.tj,sig,'-');
semilogy(r.tj,sig,'r-'); hold on; semilogy(r.tj,-sig,'-');
semilogy(r.tj,sum(abs(muall),1)/N,'g-');
%semilogy(r.tj,sqrt(sum(muall.^2,1),'g-');
axis tight; xlabel('$$t$$','interpreter','latex');
v = axis; v(3) = 1e-5; axis(v);
%ylabel('$$u_{tot}(x_0,t)$$','interpreter','latex');
title('(c) \quad total wave signal at $$x_0$$, and density 1-norm', 'interpreter','latex');
h=legend('positive $$u_{tot}(x_0,t)$$','negative $$u_{tot}(x_0,t)$$','$$\|\mu(\cdot,t)\|_1/N$$','location','south');
set(h,'interpreter','latex');
%hold on; plot(t0,utots(end),'k.','markersize',10);  % add to (c)
set(gcf,'paperposition',[0 0 5 3]);
print -depsc2 scatt_sig.eps


% ================== (d) error convergence at x_0, t0, vs np = 6,9,...
t0 = 5.0;    % 3.5 is in negative reflected peak
tes = 1:0.1:9; utes = nan(numel(tes),numel(run));   % reinterp to reg grid
for ru=1:numel(run), r=run{ru};    % loop over diff np runs
  signal = r.utot(jz+numel(gz)*(jx-1),:);   % utot at x0
  [jmax jmin w] = interpmat(-t0,r.dt,m);   % use order m+2 interp
  utots(ru) = dot(w(end:-1:1),signal(-jmax:-jmin));  % guessing offset and that w goes backwards in time!
  [jmax jmin w] = interpmat(-tes,r.dt,m);
  utes(:,ru) =  w(:,end:-1:1) * signal(-jmax:-jmin)';
  hs(ru) = r.h;
end
%errs = abs(utots(1:end-1)-utots(end));   % at single time t0
%errs = max(abs(utes(:,1:end-1)-utes(:,end)));
errs = median(abs(utes(:,1:end-1)-utes(:,end)));   % median along time axis
figure(2); clf; hold off
hs =hs(1:end-1);  sc = 0.35;
loglog(hs,errs,'+-'); hold on; plot(hs,(hs/sc).^o.p,'r--');
axis tight
xlabel('$h$','interpreter','latex');
%ylabel('empirical error in $u(x_0,t_0)$','interpreter','latex');
title('(d) \quad convergence at target $$x_0$$','interpreter','latex');
h = legend('median error in $$u_{tot}(x_0,\cdot)$$','order $$p=6$$','location','southeast');
set(h,'interpreter','latex');
%text(0.05,1e-3,'replace with better nr data, or dt=0.8 h');
set(gcf,'paperposition',[0 0 5 3]);
print -depsc2 scatt_ptconv.eps

%memorygraph('plot',bytes,est_times,cpu_times,cpu_usages,labelstrings,labeltimes);
