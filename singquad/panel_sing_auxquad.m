function [z w] = panel_sing_auxquad(targs,o)
% PANEL_SING_AUXQUAD  Auxiliary singular quadrature nodes in parameter plane
%
% [z w] = panel_sing_auxquad(targs,par) spits out auxiliary quadrature nodes and
%   weights for singular kernels with singularity no worse than 1/r on a 3x3
%   near panel neighborhood in the parameter plane
%
% Called with no arguments does a self-test
%
% Quadrilateral patches with 3x3 neighbor connectivity, fixed # aux nodes
% (q) per target in the panel, for now.
%
% Input:
%   targs (2*ntarg) target points as column vectors in parameter plane [-1,1]^2
%   opts (struct) controls options such as:
%        opts.nt, opts.nr  - # theta nodes and # radial nodes
% Output:
%   z (2 by q by ntarg) aux quadr nodes in 3x3 patch parameter plane [-3,3]^2
%   w (q by ntarg) aux quadr weights for parameter plane
%
% See also: GAUSS

% Barnett 7/11/15 based on try_local_panel_quadr.m
% todo: decide if output ordering is bad for RAM access

if nargin==0, test_panel_sing_auxquad; return; end
if nargin<2, o=[]; end
if isfield(o,'nt'), nt = o.nt; else nt = 24; end  % convergence params (# nodes)
if isfield(o,'nr'), nr = o.nr; else nr = 12; end
q = 4*nt*nr;         % total # nodes per target
[d N] = size(targs); if d~=2, error('targs must be 2*something'); end  % ntarg
z = zeros(2,q,N); w = zeros(q,N);                  % allocate outputs
[xr wr] = gauss(nr); xr = (xr(:)'+1)/2; wr = wr/2; % on [0,1], both row vecs
[xt wt] = gauss(nt); xt = (xt(:)'+1)/2; wt = wt/2; % "
for j=1:N, x = targs(:,j);       % ..... loop over targets x in R2
  corners = 3*[1 1;-1 1];  % for 3x3 panels, each col is a corner of triangle  
  indoff = 0;              % node index offset
  for tri=1:4              % triangle is CCW (x,corner(:,1),corner(:,2))
    ray1 = corners(:,1)-x; edge = corners(:,2)-corners(:,1);      
    perp = ray1 - (dot(ray1,edge)/dot(edge,edge))*edge; % perp vec from x
    t1 = atan2(corners(2,1)-x(2),corners(1,1)-x(1)); % t quadrature in theta
    t2 = atan2(corners(2,2)-x(2),corners(1,2)-x(1));
    if t2<t1, t2=t2+2*pi; end         % ensure CCW
    ts = t1 + (t2-t1)*xt;             % angles
    tperp = atan2(perp(2),perp(1));   % angle of perpendicular to edge
    rs = norm(perp) * sec(ts-tperp);  % length of rays
    rays = repmat(rs,[2 1]) .* [cos(ts);sin(ts)];  % vecs targ to edge pts
    wrays = (t2-t1) * rs.^2 .* wt;    % ray weight factors
    for i=1:nt      % loop over theta nodes
      inds = indoff + (i-1)*nr+(1:nr); % indices in full aux node list
      z(:,inds,j) = bsxfun(@plus, x, rays(:,i)*xr);     % outer prod
      w(inds,j) = wrays(i)*(xr.*wr); % beta dbeta part of polar metric
    end
    corners = [0 -1;1 0]*corners;     % rotate corners around outer square
    indoff = indoff + nt*nr;
  end  
end                              % .....
%%%%%%

function test_panel_sing_auxquad   % self-test: convergence only for now

if 1 % plot default aux nodes
  x = [.96;.9];
  [z w] = panel_sing_auxquad(x);
  figure; plot(z(1,:),z(2,:),'.'); hold on; plot(x(1),x(2),'ro');
  xp = [1;1]*[-3:2:3]; yp = [3;-3]*ones(1,4); plot([xp yp],[yp xp],'g-');
  axis equal tight
end

% targs = [0 .96;0 .9];  % simple couple of targ pts
p = 8; [x1 x2] = meshgrid(gauss(p)); targs = [x1(:)';x2(:)']; % targ panel quad
aspects = [1 2 4 8];  % toy aspect ratios in parametrization of surface
ns = 6:2:32;       % convergence in # theta nodes per triangle patch
o.verb = 2;         % verbosity: 0 just text, 1 just conv plot, 2 nodes plot
N = size(targs,2);
figure;
for n=1:numel(aspects), aspect = aspects(n)  % ======== aspect ratio sweep
  % dense and ker funcs need to vectorize over cols...
  dens = @(x) sin(1.3*x(1,:) + 0.8*x(2,:) + 0.7); % a density (.5 lambda/panel)
  % Lap SLP after simple aspect-ratio change...
  ker = @(x,y) (1/4/pi)./sqrt((aspect*(x(1,:)-y(1,:))).^2+(x(2,:)-y(2,:)).^2);
  us = nan(numel(ns),N);              % convergence of potentials
  for m=1:numel(ns)        % .......... convergence loop
    o.nt = ns(m); o.nr = round(o.nt);
    [z w] = panel_sing_auxquad(targs,o);
    for j=1:N, zj = z(:,:,j);     % make 2*q for each target
      % eval ker, dens @ aux nodes, do quadr...
      us(m,j) = sum(w(:,j)'.*ker(targs(:,j),zj).*dens(zj));
    end
    fprintf('nr=%d nt=%d (%d aux nodes / targ): \tu(1) = %.16g\n',o.nr,o.nt,size(w,1),us(m,1))
  end                     % ...........
  e = max(abs(us-repmat(us(end,:),[numel(ns) 1])),[],2);       % max abs errors
  semilogy(ns,e,'+-'); axis([min(ns) max(ns) 1e-16 1]); hold on;
  i=round(numel(ns)/2); text(ns(i),5*e(i),sprintf('aspect = %.3g',aspect))
end                 % ==========
xlabel('n_{\theta} = 2n_r'); ylabel('|u-u_{conv}|');
title('test panel sing auxquad: Lap SLP on various rectangles');
