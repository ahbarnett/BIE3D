% test t-domain GRF for wave equation in exterior of a torus.
% With Hagstrom+Greengard.
% Barnett 12/15/16

clear
% surface...
so.a=1; so.b = 0.5; [s N] = create_panels('torus',so,[]); % default # pans
[x nx] = getallnodes(s);

% surf data...
%t0 = 1.0; T = @(t) exp(-(t/t0).^2/2); Tt = @(t) (-t/t0^2).*T(t); % data t-func
w0 = 2.0; T = @(t) cos(w0*t); Tt = @(t) -w0*sin(w0*t);
eps = 1e-5; t = 0.3;       % test time
fprintf('error in time deriv: %.3g\n',(T(t+eps)-T(t-eps))/(2*eps) - Tt(t))
clear t
xs = [0.9;-0.2;0.1];   % src pt for data, must be inside

% exterior GRF test...
t.x = [1.3;0.1;0.8]; t.N = 1; % far test targ pt (a fake panel), exterior
%t.x = [-.8;.1;.2]; % ... or int GRF, should get extinction (u=0)
t.t = 2.1;   % test target time
r = sqrt(sum(bsxfun(@minus,t.x,x).^2,1));   % dists from test pt to surf pts
tret = t.t - r;   % retarded times on surf rel to GRF test pt

% now eval retarded sig, tau, tau'...
[f,fn,ft] = data_ptsrc(xs,T,Tt,tret,x,nx);
retsig = -fn; rettau = f; rettaut = ft;  % ext GRF: D.u - S.un
%showsurffunc(s,retsig); title('GRF: retarded sigma'); showsurffunc(s,rettau); title('GRF: retarded tau'); showsurffunc(s,rettaut); title('GRF: retarded tau_t');

u = LapSeval_panels(t,s,retsig) + LapDeval_panels(t,s,rettau) + RetDeval_panels(t,s,rettaut)  % u(x,t) = S.sigma + (wave eqn D).tau
uex = data_ptsrc(xs,T,Tt,t.t,t.x)
fprintf('N=%d, interior GRF test at 1 pt: u err = %.3g\n', N, u-uex)
% error 1e-8
