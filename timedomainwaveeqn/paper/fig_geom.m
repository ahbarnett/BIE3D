% paper figures for geom and spatial quadr. Barnett 11/20/18
clear; run('../../bie3dsetup')

% 3D PLOT OF SHAPES, PANELS, SOME NODES.....................
for shape=1:2        % 1=torus, 2=cruller.
  so.a=1; b = 0.5;
  wc = 0.1;  % surf modulation ampl (0 for plain torus):   cruller
  wm = 3;   % # wobbles in minor, poloidal (note swapped from 2013)
  wn = 5;   % # wobbles in toroidal, major
  f = @(t,p) b + wc*cos(wm*t+wn*p);         % wobble func, then its partials...
  ft = @(t,p) -wc*wm*sin(wm*t+wn*p); fp = @(t,p) -wc*wn*sin(wm*t+wn*p);
  so.b = {f,ft,fp};     % pass in instead of b param
  if shape==1, so.b = b; end   % override w/ plain torus, otherwise cruller

  so.np=15; so.mp = round(so.np/3*2);   % spatial discr panel numbers
  o.p = 6;                    % p = spatial "order" (G-L nodes per panel side)
  [s N] = create_panels('torus',so,o);
  [x nx] = getallnodes(s);
  showsurffunc(s,0.1+0*x(1,:));
  colorbar off; colormap(gray); axis off; hold on;
  h = showpanel(s); set(h,'color',0.6*[1 1 1]);   % adds lines
  
  k = 101; t = s{k};            % good panel to show for (np,mp) = (15,10)
  plot3(t.x(1,:),t.x(2,:),t.x(3,:),'k.');
  neipan = s(t.nei);
  [neix] = getallnodes(neipan);
  plot3(neix(1,:),neix(2,:),neix(3,:),'g.');

  view(-30,35); if shape==1, str='(a)'; else str='(b)'; end
  text(-1,1.2,1,str);
  
  set(gcf,'paperposition',[0 0 6 4]);
  % fails: EPS crash evince and acroread!!!
  %if shape==1, print -depsc2 geomtorus.eps
  %else, print -depsc2 geomcruller.eps
  %end
  if shape==1, print -dpng -r600 geomtorus.png; system('convert -trim geomtorus.png eps2:geomtorus.eps');
  else, print -dpng -r600 geomcruller.png; system('convert -trim geomcruller.png eps2:geomcruller.eps');
  end
end
