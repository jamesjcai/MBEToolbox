function [aln] = readfasta2(filename,seqtype,geneticcode)
%READFASTA - Reads data from a FASTA formatted file into a MATLAB structure
%A file with a FASTA format begins with a right angle bracket (>) and a single
%line description.
%
% Syntax:  [aln] = readfasta2('filename',seqtype,geneticcode)
%
% Inputs:
%    filename    - FASTA formatted file (ASCII text file).
%    seqtype     - Type of sequences. 1 - non-coding nucleotide; 2 - coding
%                  nucleotide; 3 - protein
%    geneticcode - Genetic code. 1 - Standard; 2 - Vertebrate Mitochondrial;
%                  3 - Yeast Mitochondrial; 4 - Mold, Protozoan, Coelenterate
%                  Mitochondrial, and Mycoplasma/Spiroplasma; 5 - Invertebrate
%                  Mitochondrial; 6 - Ciliate, Dasycladacean, and Hexamita
%                  Nuclear; 0 - non-coding
%
% Outputs:
%    aln    - alignment structure
%
% Examples:
%
% >> [aln] = readfasta2('input.fas',1,0)   % for non-coding DNA/RNA
% >> [aln] = readfasta2('input.fas',2,1)   % for protein-coding DNA
% >> [aln] = readfasta2('input.fas',3,1)   % for protein
%
% See also: READPHYLIP, WRITEFASTA

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
       {'*.fasta;*.fas', 'FASTA Format Files (*.fasta, *.fas)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a FASTA file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];

end

if nargin < 3
	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	if (isempty(seqtype)||isempty(geneticcode)), aln=[]; return; end
end



disp(['Reading ',filename]);

txtc = textread(filename,'%s','delimiter','\n','whitespace','','bufsize', 40950);
% txtc = textread('seq_examples\mhc.fas','%s','delimiter','\n','whitespace','');

mt = find(cellfun('isempty',txtc));
% eliminate empty lines
txtc(mt) = [];

txt=char(txtc);
mt=find(txt(:,1)=='>');
names=txt(mt,[2:end]);
[n,m]=size(txt);


% add one line at the end of file
mt2=cat(1,mt,n+1);
S=[];
for (k=1:length(mt)),
	sblock = txt([mt2(k)+1:mt2(k+1)-1],:);
	seq = i_mat2vector(sblock);
	S=strvcat(S,upper(seq));
end

aln.seqtype = seqtype;
aln.geneticcode = geneticcode;
aln.seqnames = i_mat2cell(names);
aln.seq=S;
aln=encodealn(aln);