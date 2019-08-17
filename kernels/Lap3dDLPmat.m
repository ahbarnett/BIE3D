function [A An] = Lap3dDLPmat(t,s)
% LAP3DDLPMAT.  dense Laplace DLP Nystrom eval matrix from sources to targets
%
% [A An] = Lap3dDLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% NB: fills M*N matrices, so assumes small cases. Faster if An not needed.
%
% Without args, does math then timing test.
% Todo: * replace with C+omp MEX.

% Barnett 8/16/19

if nargin==0, test_Lap3dDLPmat; return; end

M = size(t.x,2); N = size(s.x,2);
d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
ny = (1/4/pi) * s.nx;     % apply prefactor here, cheaper (don't use s.nx now!)
ddotsn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
A =  bsxfun(@times, ddotsn ./ (sqrt(rr).*rr), s.w);  % including src quadr wei
if nargout>1                  % targ deriv wanted ... not the fastest
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  tndotsn = bsxfun(@times,t.nx(1,:)',ny(1,:)) + bsxfun(@times,t.nx(2,:)',ny(2,:)) + bsxfun(@times,t.nx(3,:)',ny(3,:));
  An = bsxfun(@times, (1./(sqrt(rr).*rr)).*(-3./rr.*ddottn.*ddotsn + tndotsn), s.w);  % dipole deriv, incl src quad wei
end

%%%%%%%%%
function test_Lap3dDLPmat
% math test: Gauss's Law from torus-like surface...


% timing test...
ns = 5e3;
nt = 1e4;
s.x = randn(3,ns);
s.nx = randn(3,ns); s.nx = bsxfun(@times,s.nx,1./sqrt(sum(s.nx.^2,1)));
%norm(s.nx(:,1))-1
t.x = randn(3,nt);
tic;
%profile clear; profile on
A = Lap3dDLPmat(t,s);
%profile off; profile viewer
t=toc;
fprintf('filled lap dipole pot mat %d*%d in %.3g s: %.3g Gels/s\n',nt,ns,t,ns*nt/t/1e9)
