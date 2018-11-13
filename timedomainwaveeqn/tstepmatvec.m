function x = tstepmatvec(muhist, predcorr, R, wpred, o)
% TSTEPMATVEC.  Apply one time-step (zero data) to density history vector.
%
% Inputs:
%  muhist = tstepmatvec(muhist, predcorr, R, wpred, o) applies a time-step to
%  muhist - density history column vector of length n*N, 
%  predcorr = -1 (implicit), 0,1,2... use predictor with this # corrector steps
%  wpred = extrapolation vector (from extrap(m) for order m+2).
%  R - surf self-eval matrix, acts on density history vec to eval td-BIE rep.
% Output:
%  muhist - updated density history vector

% extracted from tmarch.m, Barnett 11/10/18

persistent Rnow Cnow Bnow Lnow Unow pnow iter
% see: https://www.mathworks.com/help/matlab/ref/persistent.html
% use clear tstepmatvec to reset the persistent vars

if nargin<5, o=[]; end
if ~isfield(o,'shift'), o.shift=0; end

N = size(R,1);     % # surf pts
n = size(R,2)/N;   % # time hist pts in dens hist vec and R mats

if isempty(Rnow)   % one-time build all small objects
  iter = 0;
  Rnow = R(:,n:n:end); % pull out N*N current-time mat (last time index not first!)
  Cnow = diag(Rnow) + 1/2;    % col vector, note includes the 1/2
  Bnow = Rnow - diag(diag(Rnow));   % off-diag part
  [Lnow Unow pnow] = lu(speye(N)/2 + Rnow,'vector'); % direct factor on-surf sys
end

% do the matvec...
muhist = reshape(muhist,[n N]); muhist = muhist([2:n,n],:);   % shuffle back 1
muhist(end,:)=0; muhist = reshape(muhist,[n*N 1]);       % ensure munow=0
  
% RUNS 5X SLOWER THAN IN tmarch.m ------ *** WHY?
rhs = -R*muhist;    % eval hist: RHS for "now" lin sys
if predcorr>=0                             % (NB 0 steps just extrapolates mu)
  muh = reshape(muhist,[n,N]); muh=muh(n-numel(wpred):n-1,:);
  munow = (wpred * muh)';  % predictor, kicks off iter for munow. row vec len N
  for k=1:predcorr         % corrector steps, expressed via change in munow...
    dmunow = (rhs - Rnow*munow - munow/2) ./ (Cnow + o.shift);  % shifted
    munow = munow + dmunow;
  end
else    % implicit solve
  munow = Unow\(Lnow\rhs(pnow));
end
muhist(n:n:end) = munow;           % write into "now" entries
%x =muhist; % pass out
iter=iter+1;
%fprintf('iter %d\n',iter)
