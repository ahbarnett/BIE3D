% gather stabilty results for wobbly torus known pt src,
% m=6 & various p cases (6 or 8)
% Barnett 10/23/18

clear
%mats = {'stab_wtorus_p6_m6_pulse_al1.mat'};   % p=6 data
% p=8 data:
mats = {'stab_wtorus_p8_m6_pulse_al1.mat', 'stab_wtorus_p8_m6_pulse_al1_np18.mat'};
% get full set of timesteps and other expt-wide params from run 1...
load(mats{1}); l=load(mats{2});
j=find(nps==18); % add in np=18 data (separate run): errs, gros
errs(:,j,:,:)=l.errs; gros(:,j,:,:)=l.gros;
nps = nps(1:j);    % is nothing beyond np=18!
nnp = numel(nps); ndt = numel(dts);
dxs = 2*pi./(o.p * nps);   % h=dx. mean radius 1. todo: est h min, h max.
fprintf('min err over any dx, dt : %.3g\n',min(errs(:)))  % nans don't count

iCFL = 0.5;   % sets comparison line plot slope

for ipc = 1:5          % ..................
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
