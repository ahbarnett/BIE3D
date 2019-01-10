% Extract orders for torus and cruller, the longer T=50 runs from 12/20/18.
% Barnett 1/7/19

clear
fnam = {'torus/t','wobblytorus/wt'};
nams = {'torus','cruller'};
surfareas = [19.7392088021787 22.6225981826694];
ps = [4 6 8]
lt = {'+-','.-','o-'};
figure(1);
for type=1:2     % dt-conv then dx-conv
for shape=1:2   % torus then cruller
  sp = shape+(type-1)*2;
  chr = 96 + sp;  % character code for a,b,c,...
  subplot(1,4,sp); hold off; set(gca,'ColorOrderIndex',1);
  for iii = 1:numel(ps), p=ps(iii);
  m=p-2;   % how we decided it in runs
  load(sprintf('../expts/%s_T50_m%d_p%d_pulse.mat',fnam{shape},m,p));
  nnp = numel(nps); ndt = numel(dts);
  %dxs = 2*pi./(o.p * nps);   % dx old! mean radius 1. todo: est h min, h max
  Ns = o.p^2*nps.^2/1.5;    % list of # dofs, uses mp=(2/3)np.
  dxs = sqrt(surfareas(shape)./Ns);   % defn of h = dx_typ for surface
  if type==1
    j = nnp;   % for dt-conv: pick finest h
    loglog(dts,errs(:,j),lt{iii});
    dt0 = 0.8; e0 = 0.03; if shape==2, e0=0.02; end  % empirical fits
    i=get(gca,'ColorOrderIndex'); set(gca,'ColorOrderIndex',i-1);
    hold on; plot(dts,e0*(dts/dt0).^(m+2),'--');
    xlabel('$\Delta t$','interpreter','latex');
    axis([min(dts) max(dts) 5e-10 0.03]);
    title(sprintf('(%c) %s, fixed $h=$ %.2g',chr,nams{shape},dxs(j)),'interpreter','latex')
  else                 % dx-conv:
    joff=18-3*m/2; if shape==1, joff=joff+1; end  % which dt to start on to be in h-dom (must be as in fig_stability.m)
    dts(joff-1)
    errsli = []; for j=1:nnp,errsli(j)=errs(joff-j,j); end  % carve thru errs
    loglog(dxs,errsli,lt{iii});
    h0 = 0.15; e0 = 1e-4; if p==4, e0=3e-4; end     % empirical fits
    i=get(gca,'ColorOrderIndex'); set(gca,'ColorOrderIndex',i-1);
    hold on; plot(dxs,e0*(dxs/h0).^p,'--');
    xlabel('$h$','interpreter','latex');
    axis([min(dxs) max(dxs) 5e-10 2e-4]);
    title(sprintf('(%c) %s, $\\Delta t$ on $h$-dom.\\ curve',chr,nams{shape}),'interpreter','latex')
  end
  ylabel('max error in $u(x,t)$','interpreter','latex');
  if sp==4, h=legend('$p=2q+2=4$','order 4','$p=2q+2=6$','order 6','$p=2q+2=8$','order 8');
    set(h,'location','southeast','interpreter','latex');
    %v = get(h,'position'); v(1:2)=v(1:2)+[0.02,-0.02]; set(h,'position',v);
  end
  end
end
end
set(gcf,'paperposition',[0 0 12 3.8]);
print('-depsc2',sprintf('orders.eps'))
