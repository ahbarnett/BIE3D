% gather stabilty results for wobbly torus known pt src,
% p=6, m=4
% Barnett 10/19/18

clear
mats = {'stab_wtorus_p6_m4_pulse_al1_npupto12.mat',
        'stab_wtorus_p6_m4_pulse_al1_np15+.mat',
        'stab_wtorus_p6_m4_pulse_al1_np21_serial.mat',
        'stab_wtorus_p6_m4_pulse_al1_np24_serial.mat'};
% get full set of timesteps and other expt-wide params from run 1...
l = load(mats{1}); npc = l.npc; nbe = l.nbe; nnp = 7;
dts = [0.03 0.04 l.dts]; ndt = numel(dts);
errs = nan(ndt,nnp,nbe,npc); gros=errs;  % all errors, growth facs

for r = 1:numel(mats)    % load from each run
  l = load(mats{r});
  if l.dts(1)==0.05, idt = 3:ndt;    % which dts to write into
  else, idt = 1:ndt; end
  if r==1, inp = 1:3; elseif r==2, inp = 4:5; else, inp = r+3; end
  nps(inp) = l.nps;     % the np values from each run
  tnps(inp) = l.tnps;
  errs(idt,inp,:,:) = l.errs;
  gros(idt,inp,:,:) = l.gros;
end
dxs = 2*pi./(l.o.p * nps);   % h=dx. mean radius 1. todo: est h min, h max.
ibe = 1; ipc = 3;  % ipc=npc for implicit
figure;
subplot(2,1,1);
h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), log10(errs(:,:,ibe,ipc)));
view(2); axis xy; v = axis; v([1 3])=0; axis(v);
shading interp; set(h,'edgecolor',[0 0 0]);
caxis([-10 0]); colorbar
title(sprintf('wtorus $\\log_{10}$ max error: $\\beta$ = %.3g, \\#pc = %d',l.betas(ibe),l.pcs(ipc)),'interpreter','latex');
xlabel('mean dx'); ylabel('dt');
hold on; x=[0 0.2]; h=plot3(x,0.5*x,0*x,'m--','linewidth',2);
subplot(2,1,2);
h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), gros(:,:,ibe,ipc));
view(2); axis xy; v = axis; v([1 3])=0; axis(v);
shading interp; set(h,'edgecolor',[0 0 0]);
caxis([-1 1]); colorbar
title(sprintf('wtorus growth factor: $\\beta$ = %.3g, \\#pc = %d',l.betas(ibe),l.pcs(ipc)),'interpreter','latex');
xlabel('mean dx'); ylabel('dt');
hold on; x=[0 0.2]; h=plot3(x,0.5*x,0*x,'m--','linewidth',2);

%print -depsc stab_wtorus_p6_m4_beta1_implicit.eps

figure; loglog(nps,tnps,'+'); xlabel('n_p'); ylabel('run time');
hold on; plot(nps, 0.3*nps.^4,'r--');

