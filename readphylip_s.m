function [aln]=readphylip_s(filename,seqtype,geneticcode)
%READPHYLIP_S - Reads data from a sequential PHYLIP formatted file into a MATLAB structure
%
% Syntax:  [aln]=readPhylip_s('filename',seqtype,geneticcode,'format')
%
% Inputs:
%    filename    - Phylip formatted file (ASCII text file).
%    seqtype     - Type of sequences. 1 - non-coding nucleotide; 2 - coding
%                  nucleotide; 3 - protein
%    geneticcode - Genetic code. 1 - Standard; 2 - Vertebrate Mitochondrial;
%                  3 - Yeast Mitochondrial; 4 - Mold, Protozoan, Coelenterate
%                  Mitochondrial, and Mycoplasma/Spiroplasma; 5 - Invertebrate
%                  Mitochondrial; 6 - Ciliate, Dasycladacean, and Hexamita
%                  Nuclear; 0 - non-coding
%    format      - (optional) Format of the input Phylip file: 'sequential' or
%                  'interleaved'
%
% Outputs:
%    aln    - alignment structure
%
% Examples:
%
% >> [aln] = readPhylip_s('input.phy',1,0)   % for non-coding DNA/RNA
% >> [aln] = readPhylip_s('input.phy',2,1)   % for protein-coding DNA
% >> [aln] = readPhylip_s('input.phy',3,1)   % for protein
%
% See also: READFASTA, WRITEPHYLIP

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.phylip;*.phy', 'Sequential Phylip Files (*.phylip, *.phy)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Phylip file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];
end

if nargin < 3
	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	if (isempty(seqtype)|isempty(geneticcode)), aln=[]; return; end
end


	file = fopen(filename, 'r');
	disp(['Reading ',filename]);
	[nm,x] = fscanf(file,'%d',[1 2]);
	if (x~=2) error('NOT PHYLIP FORMAT'); end
	n=nm(1);
	m=nm(2);

	% names=32*ones(n,10);	% 32 = space
	% seqs=45*ones(n,m);	% 45 = -
	temp=[];

	tline = fgetl(file);
	for (k=1:n)
	    if (feof(file)) error('NOT ENOUGH SEQUENCES'); end
	    tline = fgetl(file);
	    disp(tline(:,1:10))

	    if ~ischar(tline), break, end
		if(size(tline,1)~=0 & size(tline,2)~=0)
		    temp = cat(1,temp,tline);
		end
	end

	names={};
	for (k=1:n),
      		names{k}=removeblanks(temp(k,1:10));
	end
	% names=cellstr(temp(:,1:10));

	S=temp(:,11:end);
	S=upper(removeblanks(S));
	fclose(file);


aln.seqtype = seqtype;		% 'DNA/RNA', 'Protein', '' = Not assigned
aln.geneticcode = geneticcode;	% 0 means non-coding (0..1..13), -1 = Not assigned
aln.seqnames = names;		% arry cell {}
aln.seq=S;			% S
aln=encodealn(aln);
