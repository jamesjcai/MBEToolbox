function [h]=runphylip(aln,command,viewoutput)
%RUNADDIN - Run add-in Phylip commands
%
% Syntax: runaddin(aln,'command')
%
% Inputs:
%    aln       - Alignment structure
%    command   - 'dnapars'|'dnaml'|'protpars'|'proml'
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (nargin<3), viewoutput=1; end

dirstr=chdir2where('runphylip.m');

if (ispc),
writephylip_s(aln,[dirstr,'\infile']);
else
writephylip_s(aln,[dirstr,'/infile']);
end

switch (command)
    case ('dnapars')
	gcbox=mbe_dnapars;
    case ('dnaml')
	gcbox=mbe_dnaml;
    case ('protpars')
	gcbox=mbe_protpars;
    case ('proml')
	gcbox=mbe_proml;
    otherwise
	error('WrongOptions','Wrong options')
end

handles = guihandles(gcbox);
N=size(aln.seq,1);
X=num2str([1:N]');
X=char(aln.seqnames);
set(handles.outgrno, 'String', X);


%set(handles.figure1,'WindowStyle','modal');
if (nargout>0),
      h=handles;
end