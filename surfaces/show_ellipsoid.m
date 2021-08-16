function show_ellipsoid(E,R,t,colorvec)
% show_ellipsoid(E,R,t,colorvec)
%
% Barnett 8/15/21
if nargin<4, colorvec = .5*[1 1 1]; end
b = ellipsoid(E(1),E(2),E(3));   % baseline object at the oridin, aligned
b = setupsurfquad(b,[100 50]);
b.x = t + R * b.x;    % rot then transl, b just for vis
plot3(b.x(1,:),b.x(2,:),b.x(3,:),'.','color',colorvec);
axis equal vis3d   % gonna need always
