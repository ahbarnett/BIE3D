function [t utrg rhsnrm gnrm munrm mu] = tmarch(dt,Ttot,predcorr,gdata,R,Rtrg,wpred,o)
%
% Do explicit time marching in a t-dep BIE scheme for the wave equation.
% This routine isolates the timestepping, given self-eval and targ-eval mats,
% and RHS data.
%
% [t utrg rhsnrm gnrm munrm mu] = tmarch(dt,Ttot,predcorr,gdata,R,Rtrg,wpred,o)
%
% Inputs:
%  dt - timestep
%  Ttot - total time to evolve
%  predcorr = -1 (implicit), 0,1,2... use predictor with this # corrector steps
%  wpred = extrapolation vector (from extrap(m) for order m+2).
%  gdata - handle to Dirichlet data func of form gnow = gdata(t) returning
%            col vec of g_now on all surface points, at time t.
%  R - surf self-eval matrix, acts on density history vec to eval td-BIE rep.
%  Rtrg - target eval matrix, acts on dens hist vec to eval pot u_trg(t_now,x_trg)
%  o - optional struct, with fields such as:
%                    o.verb = 0,1,2,3 (verbosity)
%                    o.shift - diagonal shift for jacobi corrector iteration
%                    o.random - if true, random muhist (garbage) for stab check
%
% Outputs:
%  tt = time grid used (all timesteps, not just stored history by the end)
%  utrg = potential at targets on time grid
%  rhsnrm = norm (max over the boundary) of RHS in pred-corr solve on time grid
%  gnrm = norm (") of BVP bdry data on time grid
%  munrm = norm (") of density on time grid
%  mu = density on time grid (RAM-expensive to request)

% * just one targ for now -> make u a rect matrix.

% Barnett 10/12/18, 11/9/18

if nargin<8, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end
if ~isfield(o,'shift'), o.shift=0; end
if ~isfield(o,'random'), o.random=0; end
jtot = ceil(Ttot/dt);   % # timesteps
t = (1:jtot)'*dt;        % time grid (col vec)
utrg = 0*t; rhsnrm = 0*t; gnrm = 0*t; munrm = 0*t;  % col vecs
N = size(R,1);   % # surf pts
n = size(R,2)/N;   % # time hist pts in dens hist vec and R mats
keepmu = nargout>5;
if keepmu, mu = nan(N,jtot); end
Rnow = R(:,n:n:end); % pull out N*N current-time mat (last time index not first!)
if predcorr>=0   % set up for t-stepping
  Cnow = diag(Rnow) + 1/2;    % col vector, note includes the 1/2
  Bnow = Rnow - diag(diag(Rnow));   % off-diag part
else            % implicit solver (we don't propose, but is predcorr->inf limit)
  [Lnow Unow pnow] = lu(speye(N)/2 + Rnow,'vector'); % direct factor on-surf sys
%rhs = randn(N,1); mu = Unow\(Lnow\rhs(pnow)); norm(Rnow*mu + mu/2 - rhs) % check direct solve via LU
end
muhist = zeros(n*N,1);  % init dens hist vec
if o.random, muhist = randn(size(muhist)); end
t0=tic; if o.verb, fprintf('\ttmarch: starting %d steps...\n', jtot), end

for j=1:jtot           % .......................... marching loop
  muhist = reshape(muhist,[n N]); muhist = muhist([2:n,n],:); % shuffle back 1
  muhist(end,:)=0; muhist = reshape(muhist,[n*N 1]); % ensure munow=0
  
  gnow = gdata(t(j));                        % Dirichlet data vec now
  rhs = gnow - R*muhist;                     % eval hist: RHS for "now" lin sys
  if predcorr>=0                             % (NB 0 steps just extrapolates mu)
    muh = reshape(muhist,[n,N]); muh=muh(n-numel(wpred):n-1,:);
    munow = (wpred * muh)';  % predictor, kicks off iter for munow. row vec len N
    for k=1:predcorr         % corrector steps, expressed via change in munow...
      dmunow = (rhs - Rnow*munow - munow/2) ./ (Cnow + o.shift);  % shifted
      munow = munow + dmunow;
      if o.verb>2
        fprintf('\t k=%d\t ||dmunow||=%g\t ||resid||=%g\n',k,norm(dmunow),norm(Rnow*munow+0.5*munow-rhs)) % ...so can track norm
      end
      %munow = (rhs - Bnow*munow) ./ Cnow;  % plain corrector
    end
  else    % implicit
    munow = Unow\(Lnow\rhs(pnow));
  end
  muhist(n:n:end) = munow;           % write into "now" entries
  utrg(j) = Rtrg*muhist;                                    % eval pot at trg
  gnrm(j)=max(abs(gnow));
  rhsnrm(j)=max(abs(rhs));   % this is not the BVP RHS, but in the pred-corr!
  munrm(j)=max(abs(munow));
  if keepmu, mu(:,j) = munow; end
  if o.verb>1 | j==jtot
    fprintf('j=%4d (t=%.3f): |mu|=%.3g  u=%.6g\t|g|=%.3g\t|rhs|=%.3g\n',j,t(j),munrm(j),utrg(j),gnrm(j),rhsnrm(j))
  end
end                    % .........................

if o.verb, fprintf('\tdone %d steps in %.3g sec: %.3g sec per t-step\n', jtot, toc(t0), toc(t0)/jtot), end
