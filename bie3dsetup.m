function bie3dsetup
% BIE3DSETUP   puts all M-files from BIE3D package in path as absolute paths

% Barnett 7/19/16, removed mpspack 12/29/16
mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath, '/kernels']);
addpath([mfilepath, '/utils']);
addpath([mfilepath, '/test']);
addpath([mfilepath, '/singquad']);
addpath([mfilepath, '/surfaces']);
addpath([mfilepath, '/timedomainwaveeqn']);
addpath([mfilepath, '/timedomainwaveeqn/quad']);
addpath([mfilepath, '/timedomainwaveeqn/timeinterp']);
addpath([mfilepath, '/doublyperiodic']);
% hack to remove mpspack from path to avoid conflicts...
if strfind(which('showfields'),'mpspack')
  disp('note: removing mpspack from path...')
  rmpath(fileparts(which('showfields')));
end
