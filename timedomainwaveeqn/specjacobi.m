function E = specjacobi(Rnow,s,neig)
% s = corrshift
% neig even
E = nan(2,neig);
nee = neig/2;
fprintf('measure extremal eigvals of sh & un-sh corr Jacobi update matrix...\n')
t4=tic;
N = size(Rnow,1);
Cnow = diag(Rnow) + 1/2;
Bnow = Rnow - diag(diag(Rnow));   % off-diag part
E(1,1:nee) = eigs(diag(1./Cnow)*Bnow,nee,'sr');
E(1,nee+1:2*nee) = eigs(diag(1./Cnow)*Bnow,nee,'lr');
E(2,1:nee) = eigs(diag(1./(Cnow+s))*(Bnow-s*speye(N)),nee,'sr');
E(2,nee+1:2*nee) = eigs(diag(1./(Cnow+s))*(Bnow-s*speye(N)),nee,'lr');
E
fprintf('\t specjacobi done in %.3g s\n',toc(t4))
