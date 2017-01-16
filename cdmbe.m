function [pw1,pw0]=cdmbe(isconfirmed)
%CDMBE - Changes current working directory to MBEToolbox directory

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $

if nargin < 1
    isconfirmed=true;
end
pw0=pwd;
pw1=fileparts(which(mfilename));
if ~strcmp(pw0,pw1)
    [selectedButton,dlgShown]=uigetpref('MBEToolbox',... % Group
           'cdmbe_ask',...                               % Preference
           'Changing Working Directory',...              % Window title
           {'Do you want to change current working directory to MBEToolbox directory?'},...
           {'always','never';'Yes','No'},...       % Values and button strings
           'ExtraOptions','Cancel',...             % Additional button
           'DefaultButton','Yes');
    switch selectedButton
        case {'always','Yes'}
            cd(pw1);
        case {'never','No','Cancel'}
    end
end