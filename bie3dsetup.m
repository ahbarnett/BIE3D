function bie3dsetup
% BIE3DSETUP   puts all M-files from BIE3D package in path as absolute paths

% Barnett 7/19/16
mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath, '/kernels']);
addpath([mfilepath, '/utils']);
addpath([mfilepath, '/test']);
addpath([mfilepath, '/singquad']);
addpath([mfilepath, '/surfaces']);
addpath([mfilepath, '/timedomainwaveeqn']);
addpath([mfilepath, '/timedomainwaveeqn/timeinterp']);
addpath([mfilepath, '/doublyperiodic']);
