% Installation file. Adds local folders to path.

fprintf('Adding MBEToolbox folders to Matlab path... ')

%addpath(genpath(fullfile(pwd,'lib')));
%addpath(fullfile(pwd));
addpath(fileparts(which(mfilename)));

fprintf('done.\n')
disp('Type "savepath" if you wish to store changes.')
% savepath;
