% basic test for MEX interface to Hagstrom's InterpMat f90 code.
% Barnett 12/19/16
% Also see: TestInterpMat.f90

clear
dt = 0.2;
tinterp = [-1:0.05:0];  % uniform set of r=21 times in [-1,0]
m = 4;   % order 2,4,6... seems fails if odd (writes garbage to some entries)
[jmax jmin umat upmat] = interpmat(tinterp,dt,m);
'j range:'
[jmin,jmax]
ts = (jmin:jmax)*dt;  % set of grid point times
r = numel(tinterp);
'umat row sums: (should be 1)'
umat*ones(size(ts'))
'upmat row 1st-moments: (should be 1)'
upmat*ts'
figure; imagesc(ts,1:r,umat); title('umat'); hold on; plot(tinterp,1:r,'*');
xlabel('t'); ylabel('jt'); colorbar

