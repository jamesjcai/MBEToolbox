function [aln]=readphylip_i(filename,seqtype,geneticcode,noise)
%READPHYLIP_I - Reads data from an interleaved PHYLIP formatted file into a MATLAB structure
%
% Syntax:  [aln]=readPhylip_i('filename',seqtype,geneticcode)
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
% >> [aln] = readphylip_i('input.phy',1,0)   % for non-coding DNA/RNA
% >> [aln] = readphylip_i('input.phy',2,1)   % for protein-coding DNA
% >> [aln] = readphylip_i('input.phy',3,1)   % for protein
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


if (nargin<4), noise=1; end
if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.phylip;*.phy', 'Interleaved Phylip Files (*.phylip, *.phy)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Phylip file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];
end

if nargin < 3
	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	if (isempty(seqtype)||isempty(geneticcode)), aln=[]; return; end
end


	file = fopen(filename, 'r');
	if (noise), disp(['Reading ',filename]); end
	[nm,x] = fscanf(file,'%d',[1 2]);
	if x~=2, error('NOT PHYLIP FORMAT'); end
	n=nm(1);
	m=nm(2);
	fclose(file);

	txt=textread(filename,'%s','delimiter','\n','whitespace','');

	% NOT WORKING: [n,m]=strread(string(txt(1)),'%d%d','delimiter',' ');

	% remove header line
	txt(1) = [];

	% remove empty lines
	while isempty(txt{1})
	    txt(1) = [];
	end

	% find first empty string in cell array, which occurs after the first
	% consensus line
	mt = find(cellfun('isempty',txt));
	% eliminate empty lines
	txt(mt) = [];

num_seq = n;
% there are mt(1)-1 sequences
% num_seq = mt(1)-1;
%	if ~(num_seq==n)
%	    num_seq
%	    n
%    	    error('INCORRECT PHYLIP FORMAT')
%	    return;
%	end

	names={};

	for s = 1:num_seq
	    string=txt{s};
	    token = []; remainder = [];
	    token = string(1:10);
	    remainder = removeblanks(string(11:length(string)));
	    names{s}=removeblanks(token);
  	    if (noise), disp(token); end

	    for r = s+num_seq:num_seq:size(txt,1)
		% make sure that there aren't sequence numbers at the end
		remainder = [remainder removeblanks(txt{r})];
	    end

	    if ~(length(remainder)==m)
		length(remainder)
		m
		error('INCORRECT PHYLIP FORMAT')
		return;
	     end
	     S(s,:)=upper(remainder);
	end



aln.seqtype = seqtype;		% 'DNA/RNA', 'Protein', '' = Not assigned
aln.geneticcode = geneticcode;	% 0 means non-coding (0..1..13), -1 = Not assigned
aln.seqnames = names;		% arry cell {}
aln.seq=S;			% S
aln=encodealn(aln);