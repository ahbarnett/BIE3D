% 2D warm-up exercise for radial singular panel quadrature, convergence
% -------  3     3x3 square local source panels, single target x for now.
% |\   /|        Use real 2-vectors not complex for R2. Four radial triangles
% | \./ |        Explore worsening wrt aspect ratio of parametrization
% | / \ |
% |/   \|        For TDIE Greengard/Hagstrom/Epstein. Barnett 7/11/15
% ------- -3
% -3    3

clear
x = 1*[.96;.9];     % bad-case target pt (as col vec) in the [-1,1]^2 panel

aspects = [1 2 4 8];  % toy aspect ratios in parametrization of surface
ns = 6:2:32;       % convergence in # theta nodes per triangle patch
o.uniftheta = 1;    % options:  0: angle param along edge, 1: in theta
o.verb = 2;         % verbosity: 0 just text, 1 just conv plot, 2 nodes plot

for n=1:numel(aspects), aspect = aspects(n)  % ======== aspect ratio sweep

  % dense and ker funcs need to vectorize over cols...
  dens = @(x) sin(1.3*x(1,:) + 0.8*x(2,:) + 0.7); % a density (.5 lambda/panel)
  % Lap SLP after simple aspect-ratio change...
  ker = @(x,y) (1/4/pi)./sqrt((aspect*(x(1,:)-y(1,:))).^2+(x(2,:)-y(2,:)).^2);

  us = nan*ns;              % convergence of potential
  for m=1:numel(ns), nt = ns(m);     % .......... convergence loop
    nr = round(0.5*nt);              % scale the # radial nodes per triangle
    [xr wr] = gauss(nr); xr = (xr(:)'+1)/2; wr = wr/2; % on [0,1], both row vecs
    [xt wt] = gauss(nt); xt = (xt(:)'+1)/2; wt = wt/2; % "
    z = []; w = [];          % to stack up the aux nodes (2*q), wei (1*q)
    corners = 3*[1 1;-1 1];  % for 3x3 panels, each col is a corner of triangle
    for tri=1:4              % triangle is CCW (x,corner(:,1),corner(:,2))
      ray1 = corners(:,1)-x; edge = corners(:,2)-corners(:,1);      
      perp = ray1 - (dot(ray1,edge)/dot(edge,edge))*edge; % perp vec from x
      if o.uniftheta         % t quadrature in theta... (better)
        t1 = atan2(corners(2,1)-x(2),corners(1,1)-x(1));
        t2 = atan2(corners(2,2)-x(2),corners(1,2)-x(1));
        if t2<t1, t2=t2+2*pi; end         % ensure CCW
        ts = t1 + (t2-t1)*xt;             % angles
        tperp = atan2(perp(2),perp(1));   % angle of perpendicular to edge
        rs = norm(perp) * sec(ts-tperp);  % length of rays
        rays = repmat(rs,[2 1]) .* [cos(ts);sin(ts)];  % vecs targ to edge pts
        wrays = (t2-t1) * rs.^2 .* wt;    % ray weight factors
      else                   % ...or t quadrature along edge (worse)
        rays = (corners(:,2)-corners(:,1))*xt + repmat(corners(:,1) - x,[1 nt]);
        wrays = norm(perp) * norm(edge) * wt;        % ray weight factor
      end
      ztri = zeros(2,nr*nt); wtri = zeros(1,nr*nt);  % alloc for this tri
      for i=1:nt
        inds = (i-1)*nr+(1:nr);         % indices in this triangle aux node list
        ztri(:,inds) = repmat(x,[1 nr]) + rays(:,i)*xr;     % outer prod
        wtri(inds) = wrays(i)*(xr.*wr); % beta dbeta part of polar metric
      end
      z = [z ztri]; w = [w wtri];       % stack as big set of cols in R2
      corners = [0 -1;1 0]*corners;     % rotate corners around outer square
    end
    us(m) = sum(w.*ker(x,z).*dens(z)); % eval ker, dens @ aux nodes, do quadr
    fprintf('nr=%d nt=%d (%d aux nodes / targ): \tu = %.16g\n',nr,nt,numel(w),us(m))
    
    if o.verb>1 & n==1 & nt==20
      figure(4); plot(z(1,:),z(2,:),'.'); hold on; plot(x(1),x(2),'ro');
      xp = [1;1]*[-3:2:3]; yp = [3;-3]*ones(1,4); plot([xp yp],[yp xp],'g-');
      title(sprintf('aux nodes: nr=%d nt=%d uniftheta=%d',nr,nt,o.uniftheta));
      axis equal tight; hold off; end
  end               % ...............
  e = abs(us-us(end));       % error
  if o.verb, figure(3); semilogy(ns,e,'+-'); axis([min(ns) max(ns) 1e-16 1])
    hold on; xlabel('n_{\theta} = 2n_r'); ylabel('|u-u_{conv}|');
    i=round(numel(ns)/2); text(ns(i),5*e(i),sprintf('aspect = %.3g',aspect))
  end
end           % =========
if o.verb, title(sprintf('Lap SLP, rect panels, uniftheta=%d',o.uniftheta));
  hold off; end
