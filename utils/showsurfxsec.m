function h=showsurfxsec(shape,so)
% SHOWSURFXSEC  plot cross section of surface through xz (y=0) plane.
%
% h=showsurfxsec(shape,so) where shape and so is the surface params objects
%  taken by create_panels, plots the two disconnected pieces of the torus
%  slice on the current figure. Returns handles to line objects drawn.

% Barnett 1/13/19

h = [];
if strcmp(shape,'torus')
  n=200; th=2*pi*(0:n)/n;
  for phi = [0 pi]
    x = torusparam(so.a,so.b,phi*ones(size(th)),th);
    hold on; h = [h, plot(x(1,:),x(3,:),'k-')];  % pull out x-z coords
  end
else
  error('unknown shape!');
end
