function A = LapDLPmat(x,y,ny)
% Dumb dense Laplace DLP matrix fill from sources to targets
% Inputs:
% x targ (3-by-M), y src, ny src normals (both 3-by-N)
%
% NB: fills M*N matrix, so assumes small cases
%
% Without args, does timing test
% Barnett 8/1/19

if nargin==0, test_LapDLPmat; return; end

d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
rr = d1.^2+d2.^2+d3.^2;   % M*N
ny = (1/4/pi) * ny;       % apply prefactor here, cheaper
ddotn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
A =  ddotn ./ (sqrt(rr).*rr);  % kernel matrix w/o 1/4pi or src quad weights
%A = (1/4/pi) * A;

%%%%%%%%%
function test_LapDLPmat
ns = 5e3;
nt = 1e4;
y = randn(3,ns);
ny = randn(3,ns); ny = bsxfun(@times,ny,1./sqrt(sum(ny.^2,1)));
%norm(ny(:,1))-1
x = randn(3,nt);
tic;
%profile clear; profile on
A = LapDLPmat(x,y,ny);
%profile off; profile viewer
t=toc;
fprintf('filled lap dipole pot mat %d*%d in %.3g s: %.3g Gels/s\n',nt,ns,t,ns*nt/t/1e9)
