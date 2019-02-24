% Movie from fig_scattBVP.m which loads muall and utot slice from .mat file
% made by gen_scattBVPconv.m
% Barnett 1/15/19

clear
load('../expts/wobblytorus/scattBVPconv_dt0.9h.mat'); % by:gen_scattBVPconv.m
verb = 1;   % 0 for screen, 1 to write out
ru = 7; r=run{ru};   % pick run #, eg 1..7,  controls resolution np
nst = numel(r.tj);
% regen the domain
so.np = nps(ru); so.mp=round(so.np*2/3);  % panel # in major,minor directions
shape='torus'; [s N] = create_panels(shape,so,o);     % surf (torus class)
% regen the slice
dx=0.05; gx=-2:dx:2; gz = -1:dx:2; [xx zz] = meshgrid(gx,gz);
nam=sprintf('cruller_scatt_%s_pulse_m%d_p%d_np%d_hi',incwave,m,o.p,so.np);
if verb, wO=VideoWriter([nam '.avi']); wO.FrameRate=20; wO.open; end
figure; %set(gcf,'position',[300 300 800 600]);
set(gcf,'position',[100 100 1280 1024]);
nframes = nst+125;   % rotate the frozen last time step for a while
for j=1:nframes
  i = min(j,nst);    % time index
  oo=[]; oo.nofig=1; [h0 h4]=showsurffunc(s,r.muall(:,i),oo); hold on;
  h1=surf(xx,0*xx,zz,mask.*reshape(r.utot(:,i),size(zz)));
  %set(h1,'FaceLighting','none','linestyle','none','facecolor','flat');
  set(h1,'FaceLighting','none','linestyle','none','facecolor','interp');
  oo=[]; oo.dims=3; h2=showsurfxsec(shape,so,oo);   % add intersection curves
  % make caxis zoom to see the later small decaying waves...
  caxis(1.0*[-1 1]); if r.tj(i)>4, caxis(exp(0.8*(4-r.tj(i)))*[-1 1]); end
  view(-40+j/3,35); lightangle(45,0); axis tight;
  ax=gca; ax.Clipping='off'; ax.CameraPosition = ax.CameraPosition*0.6;
  %h3=title(sprintf('slice of u_{tot} and \\mu on surface: t=%.2f',tj(i)));
  h3=title(sprintf('acoustic Dirichlet scattering: total wave slice & surface density:   t=%.2f',tj(i)));
  h3.Position = [0,0,2.5];
  pos = h4.Position; h4.Position = pos + [.07,0,0,0];  % push colorbar right
  hold off; drawnow;
  if verb, writeVideo(wO,getframe(gcf)); end
  if j<nframes, clf; end   % otherwise surf not cleared, unsure why
end
if verb, close(wO);     % writes AVI movie out; now encode small MP4...
  system(sprintf('unset LD_LIBRARY_PATH; ffmpeg -i %s.avi -y -c:v libx264 -crf 20 %s.mp4',nam,nam));
end
