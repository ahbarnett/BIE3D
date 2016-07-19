function pan = add_panels_auxquad(pan, o)
% pan = add_panels_auxquad(pan, o) adds auxiliary quadr nodes for each target
%  in each of a set of panels.
%
% pan = cell array of panels.
% o = opts struct to panel_sing_auxquad
%
% Adds fields to each panel pan{k}:
%   naux - number of aux nodes
%   auxnodes (3 by naux by ntarg) - nodes in R3
%   auxnormals (3 by naux by ntarg) - node normals in R3
%   auxwei (naux by targ) - weights
%
% Notes: uses ~1000N storage

% Barnett 7/15/16
if nargin==0, test_add_panels_auxquad; return; end
[z w] = panel_sing_auxquad(pan{1}.t,o);   % same for all pans; do once (slow)
N = pan{1}.N;    % # smooth nodes (targets), same for all pans, for now
naux = size(w,1);
z = reshape(z,[2 naux*N]);  % stack all auxnodes as 2-cmpt cols
for k=1:numel(pan)
  [y ny sp] = pan{k}.chart(z);       % care: chart assumed valid in [-3,3]^2
  pan{k}.naux = naux;
  pan{k}.auxwei = pan{k}.spfac * w .* reshape(sp,[naux N]);
  pan{k}.auxnodes = reshape(y,[3 naux N]);
  pan{k}.auxnormals = reshape(ny,[3 naux N]);
end

%%%%%%%%%%%%%%%%%%%%%%
function test_add_panels_auxquad
so.a=1; so.b = 0.5; so.np = 8; so.mp = 12;
pan = create_panels('torus',so,[]);

o.nt = 16; o.nr = 8;   % # aux nodes in theta,rho
%profile clear; profile on;   % when noticed getting auxquad z,w slow
pan = add_panels_auxquad(pan,o);
%profile off; profile viewer;

W = 0;          % will be total smooth quadr wei (ie area of near panels)
figure; k=61;   % which panel to test, and plot
for knei=pan{k}.nei'    % show near pans. must be row vec for loop!
  x = pan{knei}.x; plot3(x(1,:),x(2,:),x(3,:),'g.','markersize',5); hold on;
  W = W + sum(pan{knei}.w);
end
axis equal vis3d
x = pan{k}.x; plot3(x(1,:),x(2,:),x(3,:),'r.','markersize',10);
W = W + sum(pan{k}.w);    % don't forget self for checking area
i = 21;       % which target in panel to show w/ its auxnodes
plot3(x(1,i),x(2,i),x(3,i),'k.','markersize',20);  % targ
y = squeeze(pan{k}.auxnodes(:,:,i));
plot3(y(1,:),y(2,:),y(3,:),'b.','markersize',1);  % aux
w = pan{k}.auxwei(:,i);
fprintf('diff btw nr panels area via smooth vs via aux quad = %.3g\n',sum(w)-W)
%whos
%keyboard
