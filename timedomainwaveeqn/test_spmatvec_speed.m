% some speed test for sparse matvec... simulate a bunch of targ panels

clear; p=6; N = 96*p^2;  % case p=6 on std torus
T = 1e3;  % # nodes, # time hist steps
nnzr = 8;  % nonzero els per row, ~ m (t-interp order)
p2=64;  % p^2, nodes per panel
s = randn(N*T,1);  % real dens history vec
ntargpan = 96;     % simulated (96 maxes beth RAM for p=8)
m = ntargpan*p2*N;  % # rows
i = kron(1:m,ones(1,nnzr));  % row vec of nonzero row inds
j = kron(randi(N*T,1,m),ones(1,nnzr)) + kron(ones(1,m),1:nnzr);   % contiguous col els, as in t-interp
j = mod(j,N*T)+1;   % make sure no indices fall off

% optional:
%j = randi(N*T,1,numel(i));  % non-contiguous RAM access - only 10-20% slower.

% optional:
i = mod(i-1,N)+1; m=N;  % do whole timestep as dense mat. matvec pred 0.7s all

S = sparse(i,j,randn(size(i)),m,N*T);
sparsity = numel(i)/(m*N*T)
for n=1:3, tic; t = S*s; toc, end  % repeat a couple times

% about 0.08s per matvec, per targ pan, or 0.03s per pan when 30 grouped together. 1 core (not multithreaded)
%whos
%tic; n=10; for i=1:n, t = S*s; end, toc/n % about 0.08s per matvec, per targ pan
% -> about 3s for all targ panels (would need about 20 GB RAM)

% MATLAB sp-matvec not multithreaded, but is probably RAM-limited anyway:
%http://stackoverflow.com/questions/25329654/multithreaded-sparse-matrix-multiplication-in-matlab

% conclusion: if do whole eval from NT (dens history) to N (targs) as single
% sp matvec, will take 0.7s for each of S,D, per timestep (for N=6e3, ie p=8).
% Will take around 16 GB per S,D, ie, will need desktop to run.
% Apart from RAM, good for speed!

% for p=6 panel order, still get 7 digits in GRF, and can do whole matvec
% in 0.3 s, so about 0.6 s per timestep, good!  10 GB RAM needed.



% time the shuffle of history...
clear; p=6; N = 96*p^2;  % case p=6 on std torus
T = 1e3;  % # nodes, # time hist steps
s = randn(N*T,1);
% shuffle backwards in time, leaving current time values unchanged
tic, for n=1:10, s = reshape(s,T,N); s = s([1 1:T-1],:); s = s(:); end, toc/n
% 0.007 s, for N=3400

