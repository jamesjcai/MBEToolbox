function [aln] = alignseqfile(seqtype,geneticcode,filename)
%ALIGNSEQFILE - Align sequence file.
%
% Syntax: [aln] = alignseqfile(seqtype,geneticcode,filename)
%
% Inputs:
%    seqtype       - 'DNA'|'CDS'|'PROTEIN'
%    filename      - Input sequence file (FASTA format accepted).
%    geneticcode   - Genetic code
%
% Outputs:
%    aln     - Alignment structure
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



aln=[];

if(nargin<2)
   [seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
   if (isempty(seqtype)|isempty(geneticcode)), aln=[]; return; end
end

if (nargin<3),
    [filename, pathname, filterindex] = uigetfile( ...
       {'*.*', 'Sequence Files (*.*)'}, ...
        'Select a sequence file (FASTA)');
	if ~(filename), return; end
	filename=[pathname,filename];
end



	switch (seqtype)
	    case 1
		aln = runclustalw(filename,seqtype,0);

	    case 2
		seqtype = 2;

		AlnNT = readfasta(filename,2,geneticcode);
		AlnAA = translatealn(AlnNT);
		%if (i_includeStopCodon(AlnAA))
		%    error('ERROR: Stop codon within sequences.')
		%end
		old_dir=pwd;
		cmd = 'runclustalw.m';
		%dirstr=chdir2where(cmd);
	        cd(fileparts(which(mfilename)));
		cd('addins');
	        cd('clustalw')
		i_writeFASTA(AlnAA,'infile');
		[AlnAA] = runclustalw('infile',3,geneticcode);
		aln=i_alignCDS2Protein(AlnNT,AlnAA);
		cd(old_dir);
	   case 3
		aln = runclustalw(filename,seqtype,geneticcode);
	   otherwise
		   error('Wrong seqtype.')
	   return;
	end



%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [y] = i_includeStopCodon(AlnAA)
y=0;
[n,m]=size(AlnAA.seq);

[NT,AA] = seqcode;
% AA = 'ARNDCQEGHILKMFPSTWYV*-';
removeable =find(AA=='-'|AA=='*');
Seq = AlnAA.seq;
for i=1:n,
   s=Seq(i,:);
for (j=length(s):-1:1),
	if (ismember(s(j),removeable)),
	  s(j)=[];
	else
	  break;
	end
end
	x=find(AA=='*');
	if (any(s==x)),
		y=1;
		return;
	end
end



%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function i_writeFASTA(aln,filename)
if nargin < 2
	error('No filename')
end

[n,m]=size(aln.seq);
[NT,AA] = seqcode;

if ~(aln.seqtype==3)
	error('Must be protein sequences')
else
	Seq = AA(aln.seq);
end

fid = fopen(filename,'wt');

if fid == -1,
	error('Unable to open file.')
end

for i=1:n,
   name=char(aln.seqnames(i));
   s=Seq(i,:);
for (j=length(s):-1:1),
	if (strcmp(s(1,j),'-'));
	  s(j)=[];
	else
	  break;
	end
end
   fprintf(fid, ['>%s\n'],name);
   fprintf(fid, ['%s\n'],char(s));
end
fclose(fid);


%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [aln] = i_alignCDS2Protein(AlnNT,AlnAA)
[n,m]=size(AlnAA.seq);
aln=copyalnheader(AlnNT);
S=ones(n,m*3)*5;

SA = AlnAA.seq;
SN = AlnNT.seq;

aagap=i_getcode4gap('PROTEIN');

for (i=1:n),
x=1;
y=1;
for (j=1:m),
	if ~(SA(i,j)==aagap)
		S(i,x)=SN(i,y);
		S(i,x+1)=SN(i,y+1);
		S(i,x+2)=SN(i,y+2);
		x=x+3;
		y=y+3;
	else
		x=x+3;
	end
end
end
aln.seq = S;