% fig for paper: 2D PLOT OF AUX QUAD SCHEME on std panel. Barnett 11/22/18.
clear; run('../../bie3dsetup')

p=6;
[x w] = gauss(p);   % Gauss-Legendre nodes & weights on [-1,1]
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
ww = w(:)*w; ww = ww(:)';                     % 1*p^2 weights
figure; plot(xx(1,:),xx(2,:),'.','markersize',10,'color',0.4*[1 1 1]);
hold on;
xp = [1;1]*[-3:2:3]; yp = [3;-3]*ones(1,4);  % panel grid
plot([xp yp],[yp xp],'-','color',0.7*[1 1 1]);
axis equal tight; xlabel u; ylabel v;
j = 3;   % which node in p*p grid, to show aux quad for
o.nr = 10; o.nt = 2*o.nr;     % aux quad orders
[z w] = panel_sing_auxquad(xx(:,j),o);

c = [1,.5,.5]; plot([-3 xx(1,j) -3],[3 xx(2,j) -3],'-','color',c); % triangles
plot([3 xx(1,j) 3],[3 xx(2,j) -3],'-','color',c);
plot([3 3 -3 -3 3],[3 -3 -3 3 3],'-','color',c);

plot(z(1,:),z(2,:),'.','markersize',5,'color',[0.7 0 0]);   % aux nodes

plot(xx(1,j),xx(2,j),'k.','markersize',10);
text(-4,3,'(c)','fontsize',16);
set(gcf,'paperposition',[0 0 3.5 3]);
print -depsc2 auxnodes.eps
system('cp auxnodes.eps ~/physics/leslie/hagstrom/td-acoustics/')
