% version of fig_scattBVP.m which loads muall and utot slice from .mat file
% made by gen_scattBVPconv.m
% and interpolates muall and utot to arbitrary times.
% Barnett 1/15/19

clear
%load('../expts/wobblytorus/scattBVPconv_dt0.9h.mat'); % by:gen_scattBVPconv.m
load('../expts/wobblytorus/scattBVPconv_wider_dt1.0h.mat'); % by:gen_scattBVPconv.m
ru = 5; r=run{ru}; % pick run #, eg 1..7
% regen the domain
so.np = nps(ru); so.mp=round(so.np*2/3);  % panel # in major,minor directions
shape='torus'; [s N] = create_panels(shape,so,o);     % surf (torus class)
% slice must match...
dx=0.05; gx=-2:dx:2; gz = -1:dx:2; [xx zz] = meshgrid(gx,gz);
t.x = [xx(:)';0*xx(:)';zz(:)']; t.N=numel(xx);           % x-z plane of targs
jx=37; jz=41; j=jz+numel(gz)*(jx-1); x0=xx(j); z0=zz(j);  % spatial pt
fprintf('test pt x0=(%g,%g,%g)\n',x0,0,z0);

t0s = [3.0, 7.0];           % two times, for plots (a) and (b)
for ii=1:2
  figure; t0=t0s(ii);
  [jmax jmin w] = interpmat(-t0,r.dt,m);   % use order m+2 interp
  inds = -jmax:-jmin;  w = w(end:-1:1);     % interp backwards in time!
  muall0 = r.muall(:,inds)*w(:);     % do the interp to each row of muall
  utot0 = r.utot(:,inds)*w(:);
  oo=[]; oo.nofig=1; [h0 h4]=showsurffunc(s,muall0,oo); hold on;
  h1=surf(xx-dx/2,0*xx,zz-dx/2,mask.*reshape(utot0,size(zz))); % 1/2-pix shift!
  set(h1,'FaceLighting','none','linestyle','none','facecolor','flat');
  oo=[]; oo.dims=3; h2=showsurfxsec(shape,so,oo);   % add intersection curves
  plot3(x0,0,z0,'.','markersize',20);
  caxis(1.0*[-1 1]); if ii==2, caxis(max(abs(utot0))*[-1 1]); end
  view(20,35); lightangle(45,0); axis tight;
  ax=gca; ax.Clipping='off'; ax.CameraPosition = ax.CameraPosition*0.65;
  h3=title(sprintf('(%c)     u_{tot} on \\{y=0\\} and \\mu on \\Gamma     t=%.g',ii+96,t0));
  h3.Position = [0,0,2.5];
  pos = h4.Position; h4.Position = pos + [.07,0,0,0];  % push colorbar right
  set(gcf,'paperposition',[0 0 5.8 5]);
  nam = sprintf('scatt_%d',ii);
  print('-dpng','-r300', [nam '.png']);
  system(['convert -trim ' nam '.png eps2:' nam '.eps']);
end
