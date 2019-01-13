% print some stability run total data stats.
% Barnett 1/9/19

clear
fnam = {'torus/t','wobblytorus/wt'};
ps = [4 6 8];
runs = 0; maxN = 0; minN = inf; cputot=0;
for shape=1:2
  for iii = 1:numel(ps), p=ps(iii);
  m=p-2;   % how we decided to runs
  load(sprintf('../expts/%s_T50_m%d_p%d_pulse.mat',fnam{shape},m,p));
  nnp = numel(nps); ndt = numel(dts);
  %dxs = 2*pi./(o.p * nps);   % dx old! mean radius 1. todo: est h min, h max
  Ns = o.p^2*nps.^2/1.5;    % list of # dofs, uses mp=(2/3)np.
  runs = runs + nnp*ndt;
  %max(nps)
  if max(Ns)>maxN, maxN = max(Ns); end
  if min(Ns)<minN, minN = min(Ns); end
  cputot=cputot + sum(tnps);
  max(bytes)/1e9
  end
end
runs
maxN
minN
cputot
cputot/3600
cputot/3600/24
