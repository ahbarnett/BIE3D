function I = surfinterpmat(sf,si)
% SURFINTERPMAT  Spectral interpolation matrix btw two global surface quads, 3D.
%
% I = surfinterpmat(sf,si) returns sf.N by si.N matrix interpolating smooth
%  functions sampled on the nodes on initial surface si to those on the final
%  surface sf. The surfaces must be of the same type (topology, ie, torus or
%  sphere), with grids described by s.topo, s.Nv, s.Nu, and (if sphere type)
%  s.v. The 3D locations of nodes, or surface shape, is *not* needed nor used;
%  It is mapping purely between two copies of (u,v) parameter space.
%
%  Note: is a wrapper to appropriate interpmat routines.

% Barnett 9/5/19. Fixed to m-odd sphere interp variant, 12/11/19
if iscell(si) || iscell(sf), error('si & sf cannot be patch array surfaces!'); end
if sf.topo~=si.topo, error('si and sf must be of same topology!'); end
if si.topo=='t'      % torus type
  I = peri2dspecinterpmat([sf.Nu,sf.Nv],[si.Nu,si.Nv]);
elseif si.topo=='s'  % sphere type
  if numel(si.Nu)==1 || numel(sf.Nu)==1
    error('tensor-prod spheres not implemented');
  end
  I = intxperiinterpmat(sf.Nu,sf.v,si.Nu,si.v,1);   % Nu both are vectors here
end
