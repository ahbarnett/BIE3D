function [x nx w] = getallnodes(pan)
% GETALLNODES   concatenate horizontally all nodes, normal, etc, in dof order
%
% [x nx w] = getallnodes(pan) where pan is cell array of panels returns
%  panel (ie dof) ordered nodes x, normals nx, and weights w, stacked
%  horizontally as 3-by-N or 1-by-N arrays.
%
% See:
% http://stackoverflow.com/questions/14882704/how-can-i-access-all-field-elements-of-a-structure-array-nested-in-a-cell-array

% Barnett 7/19/16
if ~iscell(pan), pan = {pan}; end
x = cellfun(@(p)p.x,pan,'uniformoutput',0); x = horzcat(x{:});
if nargout>1
  nx = cellfun(@(p)p.nx,pan,'uniformoutput',0); nx = horzcat(nx{:});
end
if nargout>2
  w = cellfun(@(p)p.w,pan,'uniformoutput',0); w = horzcat(w{:});
end
