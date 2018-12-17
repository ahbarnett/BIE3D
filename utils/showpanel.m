function h = showpanel(s,varargin)
% SHOWPANEL.  Plot a panel, or list of panels, onto current figure
%
% h = showpanel(s) where s is a panel struct plots the panel boundary
%  on current figure, returning plot handle to lines.
%  If s is cell array of panels, shows all, and returns col vec of handles.

% Barnett 11/20/18

if nargin==0, test_showpanel, return; end

if numel(s)>1      % handle cell array of panels
  h = [];
  for i=1:numel(s)
    h = [h; showpanel(s{i},varargin{:})];
  end
  return
end

% from now s is one panel
if nargin>1, o = varargin{1}; else o = []; end
if ~isfield(o,'linetype'), o.linetype = 'g-'; end

t=-1:0.01:1;    % each edge in a square loop around [1,1]^2...
x = s.chart([t, 1+0*t, t(end:-1:1), -1+0*t; -1+0*t, t, 1+0*t, t(end:-1:1)]);
h = plot3(x(1,:),x(2,:),x(3,:),o.linetype);

%%%%%%%%%%%%%
function test_showpanel        % show all panels on torus surf
so.a=1; so.b = 0.5;
so.np=12; so.mp = round(so.np/3*2);   % spatial discr panel numbers
o.p = 6;                       % p = spatial "order" (G-L nodes per panel side)
[s N] = create_panels('torus',so,o);
[x nx] = getallnodes(s);
showsurffunc(s,0.1+0*x(1,:)); colorbar off; colormap(gray); axis off; hold on;
h = showpanel(s); set(h,'color',0.6*[1 1 1]);
