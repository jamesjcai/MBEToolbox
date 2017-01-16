function [aln] = readclustal(filename,seqtype,geneticcode)
%READCLUSTAL - Reads data from a CLUSTAL formatted file into a MATLAB structure
%A file with a FASTA format begins with a right angle bracket (>) and a single
%line description.
%
% Syntax:  [aln] = readclustal('filename',seqtype,geneticcode)
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
% >> [aln] = readclustal('input.aln',1,0)   % for non-coding DNA/RNA
% >> [aln] = readclustal('input.aln',2,1)   % for protein-coding DNA
% >> [aln] = readclustal('input.aln',3,1)   % for protein
%
% Here is an example of a multiple alignment in CLUSTAL W format:
%
%CLUSTAL W (1.82) multiple sequence alignment
%
%
%FOSB_MOUSE      MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA 60
%FOSB_HUMAN      MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA 60
%                ************************************************************
%
%FOSB_MOUSE      ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS 120
%FOSB_HUMAN      ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPVVDPYDMPGTSYSTPGMSGYSSGGASGS 120
%                ********************************.***************:*.**:******
%
%FOSB_MOUSE      GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT 180
%FOSB_HUMAN      GGPSTSGTTSGPGPARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT 180
%                ****** ***** .**********************************************
%
%FOSB_MOUSE      DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD 240
%FOSB_HUMAN      DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD 240
%                ************************************************************
%
%FOSB_MOUSE      LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY 300
%FOSB_HUMAN      LPGSAPAKEDGFSWLLPPPPPPPLPFQTSQDAPPNLTASLFTHSEVQVLGDPFPVVNPSY 300
%                ****:.******.**************:*:**************************.***
%
%FOSB_MOUSE      TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL 338
%FOSB_HUMAN      TSSFVLTCPEVSAFAGAQRTSGSDQPSDPLNSPSLLAL 338
%                ***********************:**************
%
% See also: READPHYLIP, READFASTA

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
       {'*.clustal;*.aln', 'CLUSTAL Format Files (*.clustal, *.aln)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a CLUSTAL file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];

end

if nargin < 3
	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	if (isempty(seqtype)||isempty(geneticcode)), aln=[]; return; end
end
disp(['Reading ',filename]);

% txt = textread(filename,'%s','delimiter','\n','whitespace','');

fid=fopen(filename,'r');
t=textscan(fid,'%s','Delimiter','\n');
txt=t{1};
fclose(fid);

% check header line, confirm that it's a ClustalW file
if ~strncmpi(txt{1},'CLUSTAL',7)
    warning('Header does not match CLUSTAL format')
end

% remove header line and empty lines
txt(1) = [];
while isempty(txt{1}), txt(1) = []; end

% find first empty string in cell array, which occurs after the first
% consensus line
mt = find(cellfun('isempty',txt));

% eliminate empty lines
txt(mt) = [];

hasconsline=true;

if hasconsline
    % the first consensus line is in mt(1)-1
    % there are cons_loc-1 sequences
    cons_loc = mt(1)-1;
    num_seq = cons_loc-1;
else
    cons_loc = mt(1)-1;
    num_seq = mt(1)-1;
end

alnstr='';
names={};

for s = 1:num_seq,
    % make the name into a MATLAB-acceptable variable name
    name = cleantext(strtok(txt{s}),{'|','_';'.',''});
    names{s}=name;
    seqstr='';
    xxx='';
    for r = s:cons_loc:size(txt,1),
        idx=find(txt{r}==' ');
        idx=idx(1);
        % make sure that there aren't sequence numbers at the end
        xxx = [xxx strtrim(strtok(txt{r}(idx:end),'0123456789'))];
    end
    seqstr=strvcat(seqstr,char(xxx));
    try
    alnstr=[alnstr;seqstr];
    catch ME
        seqstr
        return;
    end
end


aln=struct;
aln.seqtype = seqtype;
aln.geneticcode = geneticcode;
aln.seqnames = names;
%aln.seq=nt2int(alnstr);
aln.seq=upper(alnstr);
aln=encodealn(aln);


