% figs for stability, torus and cruller, the longer T=50 runs from 12/20/18.
% Barnett 1/5/19
% Added h-dominated error extraction curves plot, 1/8/19

clear
fnam = {'torus/t','wobblytorus/wt'};
nams = {'torus','cruller'};
surfareas = [19.7392088021787 22.6225981826694];
iCFLs = [0.8, 0.7];
ps = [4 6 8];
for shape=1:2
  figure;
  for iii = 1:numel(ps), p=ps(iii);
  m=p-2;   % how we decided to runs
  load(sprintf('../expts/%s_T50_m%d_p%d_pulse.mat',fnam{shape},m,p));
  nnp = numel(nps); ndt = numel(dts);
  %dxs = 2*pi./(o.p * nps);   % dx old! mean radius 1. todo: est h min, h max
  Ns = o.p^2*nps.^2/1.5;    % list of # dofs, uses mp=(2/3)np.
  dxs = sqrt(surfareas(shape)./Ns);   % defn of h = dx_typ for surface
  fprintf('min err over any dx, dt : %.3g\n',min(errs(:)))  % nans don't count
  sp=subplot(1,numel(ps),iii);
  %h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), log10(errs));
  h = surf(ones(ndt,1)*dxs, dts'*ones(1,nnp), zeros(ndt,nnp), double(errs>1));
  xlabel('$h$','interpreter','latex');
  ylabel('$\Delta t$','interpreter','latex');
  set(gca,'xscale','log','yscale','log');     % is better than lin
  set(gca,'ytick',[.02 .03 .05 .07 .1 .15 .2 .3 .5 .8]);
  view(2); axis xy tight; v = axis; %v([1 3])=0; % lin only
  v(4)=0.5; %2.5*max(dxs);
  if shape==1, v(3)=0.035; end
  axis(v); grid off
  shading interp; set(h,'edgecolor',0.8*[1 1 1]);
  set(h,'faceAlpha','interp','faceAlpha','interp','AlphaData',errs>1,'AlphaDataMapping','scaled')
  %colormap(sp,goodbw(256))  %colormap(gray(256))
  colormap(jet(256));
  hold on;
  joff=18-3*m/2; if shape==1, joff=joff+1; end % which dt to start on to be in h-dom (must be as in fig_orders.m)
  plot(dxs,dts(joff-(1:nnp)),':','color',.5*[1 1 1]);    % h-dom curve.
  dxplot=[min(dxs) 0.1];       % add approx inv-CFL stability line
  plot(dxplot,iCFLs(shape)*dxplot,'r--','linewidth',2);
  [C,hc] = contour(dxs,dts,-log10(errs),0:10,'k-');
  clabel(C,hc,'color','k','labelspacing',50);
  chr = 96 + (shape-1)*3+iii;  % character code for a,b,c,...
  text(v(1)/(v(2)/v(1))^0.2,v(4)/(v(4)/v(3))^0.05,sprintf('(%c)',chr),'fontweight','bold');
  title(sprintf('%s $p=%d$, $2q=%d$ : $-\\log_{10}$(max err)',nams{shape},o.p,m),'interpreter','latex')
  end
  set(gcf,'paperposition',[0 0 12 5.9]);
  print('-depsc2',sprintf('stab_%s.eps',nams{shape}))
end
