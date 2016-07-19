function L = tensorprod_interp(t,s1,s2)
% TENSORPROD_INTERP   build interpolation matrix from 2D tensor-prod to arb pts
%
% L = tensorprod_interp(t,s1,s2)
% Inputs:    t - (2-by-m) target pts in some rectangle
%            s1,s2 - (length p vectors) nodes in each direction for interp from
% Outputs:  L - (m-by-p^2) dense interpolation matrix

% bits taken from Helsing-style interpmatrix from 2013
% Barnett 7/17/16

if nargin==0, test_tensorprod_interp; return; end
L1 = interpmat_1d(t(1,:),s1);   % interp to 1-coords
L2 = interpmat_1d(t(2,:),s2);   %       "   2
p1 = numel(s1); p2 = numel(s2);
m = size(t,2);
L = zeros(m,p1*p2);
for i=1:m
  L(i,:) = kron(L1(i,:),L2(i,:));  % gotta be a better way (& guessed order)
end

%%%%%%%%%%
function test_tensorprod_interp
x = gauss(16);
f = @(x) sin(x(1,:) + 1.5*x(2,:) + 0.7);  % smooth, maps cols in R2 to vals
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];   % 2-by-p^2
data = f(xx)';            % func sampled on smooth (src) nodes, col vec
t = 2*rand(2,10000) - 1;    % large # pts in [-1,1]^2
uex = f(t)';     % col vec
tic
L = tensorprod_interp(t,x,x);
toc
%whos
u = L * data;
fprintf('max abs err for interp in [-1,1]^2 : %.3g\n',max(abs(u - uex)))
