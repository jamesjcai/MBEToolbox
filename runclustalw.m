function [aln] = runclustalw(filename,seqtype,geneticcode)
%RUNCLUSTALW - Mulitple sequence alignment.
%
% Syntax: [aln] = clustalw(filename,seqtype,geneticcode)
%
% Inputs:
%    filename      - Input sequence file (FASTA format accepted).
%    seqtype       - 1 or 2 = 'DNA'; 3 = 'PROTEIN'
%    geneticcode   - Genetic code
%
% Outputs:
%    aln      - Alignment structure
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


if (nargin<3), geneticcode=1; end
if (nargin<2),
   error('MissingInput','Please input SEQTYPE.')
else
	switch (seqtype)
	    case (1)
		 stype = 'DNA';
	    case (2)
		 stype = 'DNA';
	    case (3)
		 stype = 'PROTEIN';
	   otherwise
		error('WrongSeqtype','Wrong SEQTYPE.')
	end
end

if (ispc),
	cmd = 'clustalw.exe';
	sep = '\';
else
	cmd = 'clustalw';
	sep = '/';
end

if ~(exist(filename,'file')==2),

	error('MissingInput','Missing input file.')
end

oldpath=pwd;
cdmbe;
cd 'addins';
cd 'clustalw';





if ~(exist(cmd,'file')==2),
	   error('MissingExe','Missing CLUSTALW.')
end

%h = waitbar(0.15,'Please wait...');

%dirstr=chdir2where('runclustalw.m');
%outfile = [dirstr,sep,'outfile'];
outfile = 'outfile';
optionstr=sprintf(['-align -infile="%s" -outfile="%s" -type=%s '],filename,outfile,stype);
optionstr=[optionstr, '-output=phylip -outorder=input -outputtree=phylip'];

%for i=15:30, waitbar(i/100,h); end

cmdstr = ['.', sep, cmd, ' ', optionstr, ' '];
disp(cmdstr)
[s,w] = system(cmdstr);

%for i=30:100, waitbar(i/100,h); end

if ~(s==0)
    disp(cmd)
    error('RuningTimeError','Error occurred when running CLUSTALW.');
else
    disp(w)
    aln=readphylip_i(outfile,seqtype,geneticcode);
end
cd(oldpath);

%waitbar(1,h);
%close(h);

%i_cleanfile;


function i_cleanfile()
	files={'infile.dnd','outfile','infile'};
	for (k=1:length(files)),
	      file=char(files{k});
	      if (exist(file,'file')==2),
		delete(file);
	      end
    end

