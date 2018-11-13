% gather stabilty results for torus w/ known pt src data.
% Barnett 11/09/18

clear
mats = {'tm_stab_p4_m2.mat','tm_stab_p6_m2.mat','tm_stab_p6_m4.mat','tm_stab_p6_m6.mat'};
for expt=

  *** to finish
  
% get full set of timesteps and other expt-wide params from run 1...
load(mats{1});
j=find(nps==18); % add in np=18 data (separate run): errs, gros
errs(:,j,:,:)=l.errs; gros(:,j,:,:)=l.gros;
nps = nps(1:j);    % is nothing beyond np=18!
nnp = numel(nps); ndt = numel(dts);
dxs = 2*pi./(o.p * nps);   % h=dx. mean radius 1. todo: est h min, h max.
fprintf('min err over any dx, dt : %.3g\n',min(errs(:)))  % nans don't count

iCFL = 0.5;   % sets comparison line plot slope

for ipc = 1:          % ..................
ibe = 1; %ipc = 2;  % ipc=npc for implicit
figure; set(gcf,'position',[ipc*400-400 100 450 900]);
subplot(2,1,1);
h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), log10(errs(:,1:nnp,ibe,ipc)));
view(2); axis xy; v = axis; v([1 3])=0; axis(v);
shading interp; set(h,'edgecolor',[0 0 0]);
caxis([-10 0]); colorbar
title(sprintf('wtorus $\\log_{10}$ max error: $\\beta$ = %.3g, \\#pc = %d',betas(ibe),pcs(ipc)),'interpreter','latex');
xlabel('mean dx'); ylabel('dt');
hold on; x=[0 0.2]; h=plot3(x,iCFL*x,0*x,'m--','linewidth',2);
subplot(2,1,2);
h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), gros(:,1:nnp,ibe,ipc));
view(2); axis xy; v = axis; v([1 3])=0; axis(v);
shading interp; set(h,'edgecolor',[0 0 0]);
caxis([-1 1]); colorbar
title(sprintf('wtorus growth factor: $\\beta$ = %.3g, \\#pc = %d',betas(ibe),pcs(ipc)),'interpreter','latex');
xlabel('mean dx'); ylabel('dt');
hold on; x=[0 0.2]; h=plot3(x,iCFL*x,0*x,'m--','linewidth',2);
end               % ....................

%print -depsc stab_wtorus_p6_m4_beta1_implicit.eps

if 0
  tnps = [tnps nan(1,nnp-numel(tnps))];
  figure; loglog(nps,tnps,'+'); xlabel('n_p'); ylabel('run time');
  hold on; plot(nps, 0.3*nps.^4,'r--');
end
