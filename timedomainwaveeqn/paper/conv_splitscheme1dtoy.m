% test empirical convergence order of 1D toy panel scheme w/ end singularity
% where the last panel is handled by special scheme (here just analytic!)
% The issue of nearest smooth panel formally low-order.  Barnett 1/4/19

clear
p=7;                            % Legendre pts per panel
gam = 1/2;                      % power of singular integral
% plain power law - always dominated by near sing
f = @(x) x.^-gam;               % integrand;  integrable on (0,1) for gam<1
F = @(x) x.^(1-gam)/(1-gam);    % exact indefinite integral from 0 to x
k = 50;    % add oscillatory but smooth part to see smooth conv regime...
f = @(x) x.^-gam + cos(k*x);
F = @(x) x.^(1-gam)/(1-gam) + sin(k*x)/k;

Iex = F(1);
ns = 2:2:30;

[x w] = gauss(p); x = (1+x)/2; w = w/2;   % rule on [0,1]
hs =1./ns;
Is = nan*ns;        % numerical integrals
for i=1:numel(ns), n=ns(i);
  I = F(1/n);   % "special" rule for first panel
  for j=1:n-1, I = I + (w*f((j+x)/n))/n; end
  fprintf('n=%d \t I=%.16g \t error=%.3g\n',n,I,I-Iex)
  Is(i) = I;
end
figure; loglog(hs,abs(Is-Iex),'+-'); hold on;
prefac = 0.01 * (2*p)^(-p); plot(hs, prefac*hs,'r:');  % 1st ord
axis tight; v = axis;
plot(hs, hs.^(2*p),'r--');
axis(v); xlabel('h = 1/n'); ylabel('integral error')

% indeed shows order-2p for a while, saturating at order less than 1,
% at epsilon typically 2p^-p, or maybe exp(-cp).

