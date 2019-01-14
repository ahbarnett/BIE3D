function h=showsurfxsec(shape,so,o)
% SHOWSURFXSEC  plot cross section of surface through xz (y=0) plane.
%
% h=showsurfxsec(shape,so) where shape and so is the surface params objects
%  taken by create_panels, plots the two disconnected pieces of the torus
%  slice on the current figure. Returns handles to line objects drawn.
%
% h=showsurfxsec(shape,so,o) allows options such as:
%   o.dims (default 2): if 2, plot in xy plane; if 3 do plot3 in xz (y=0).

% Barnett 1/13/19

if nargin<3, o=[]; end
if ~isfield(o,'dims'), o.dims=2; end

h = [];
if strcmp(shape,'torus')
  n=200; th=2*pi*(0:n)/n;
  for phi = [0 pi]
    x = torusparam(so.a,so.b,phi*ones(size(th)),th);
    hold on;
    if o.dims==2
      h = [h, plot(x(1,:),x(3,:),'k-')];  % pull out x-z coords in xy plane
    elseif o.dims==3
      h = [h, plot3(x(1,:),0*x(1,:),x(3,:),'k-')];  % curves in 3d space
    end
  end
else
  error('unknown shape!');
end
