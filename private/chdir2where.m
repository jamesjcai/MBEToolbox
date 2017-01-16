function [dirstr]=chdir2where(cmd)
% CHDIR2WHERE             Jump to a specified file
%
% CHDIR2WHERE filename changes MATLAB current directory to the one containing
% filename.  filename must be on the MATLAB path.  goto uses the first
% instance of filename that it finds (a la WHICH).
%
% SEE ALSO:  WHICH

% Matlab Molecular Biology & Evolution Toolbox, (C) 2004,
% Written by James J. Cai

cmd_full = which(cmd);
% 5/17/2004 Replace dirstr = cmd_full(1:size(cmd_full,2)-size(cmd,2)); with dirstr = fileparts(cmd_full);
% based on: Scott Hirsch, shirsch@mathworks.com, http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=3145&objectType=file
dirstr = fileparts(cmd_full);
cd(dirstr);