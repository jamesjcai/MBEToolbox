function [aln] = readmaf(filename,seqtype,geneticcode)
%READMAF - Reads multiple alignment format file
%
% Syntax:  [aln] = readmaf('filename',seqtype,geneticcode)
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
% Outputs:
%    aln    - alignment structure
%
% Examples:
%
% >> [aln] = readmaf('input.aln',1,0)   % for non-coding DNA/RNA
% >> [aln] = readmaf('input.aln',2,1)   % for protein-coding DNA
% >> [aln] = readmaf('input.aln',3,1)   % for protein
%
% REF: http://genome.ucsc.edu/FAQ/FAQformat#format5
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
       {'*.maf', 'MAF Format Files (*.maf)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a MAF file');
	if ~(filename), aln=[]; return; end
	filename=[pathname,filename];

end

if nargin < 3
	[seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	if (isempty(seqtype)||isempty(geneticcode)), aln=[]; return; end
end
disp(['Reading ',filename]);

txt = textread(filename,'%s','delimiter','\n','whitespace','','bufsize',40950*10);

% check header line, confirm that it's a ClustalW file
%if ~strncmpi(txt{1},'CLUSTAL',7)
%    warning('Header does not match CLUSTAL format')
%end

% remove header line
txt(1:7) = [];

% remove empty lines
while isempty(txt{1})
    txt(1) = [];
end

% find first empty string in cell array, which occurs after the first
% consensus line
mt = find(cellfun('isempty',txt));

% eliminate empty lines
txt(mt) = [];

if mod(size(txt,1),3)>0
    error('something wrong')
end


% the first consensus line is in mt(1)-1
cons_loc = mt(1)-2;

num_seq = mt(1)-2;

alnstr='';
names={};

for s = 2:num_seq:size(txt,1)
    % make the name into a MATLAB-acceptable variable name
    [temp,x]=strtok(txt{s});
    [name,x] = strtok(x);
    names{s}=cleantext(name,{'|','_';'.',''});
    seqstr='';
    atempstr='';
    for r = s:cons_loc:size(txt,1),
        idx=find(txt{r}==' ');
        idx=idx(1);
        % make sure that there aren't sequence numbers at the end
        atempstr = [atempstr strtrim(strtok(txt{r}(idx:end),'0123456789'))];
    end
    seqstr=strvcat(seqstr,char(atempstr));
    alnstr=[alnstr;seqstr];
end


aln=struct;
aln.seqtype = seqtype;
aln.geneticcode = geneticcode;
aln.seqnames = names;
aln.seq=nt2int(alnstr);
%aln.seq=upper(alnstr);
%aln=encodealn(aln);




















function txt = cleantext(txt,subs)
% CLEANTEXT converts string by removing or substituting characters
% C_TEXT = CLEANTEXT(T) replaces any characters which can't be used in
%   MATLAB variable names with an underscore, '_'.  CLEANTEXT also removes
%   all whitespace characters.  It prepends the string with 'var_', if the
%   first character left in the string isn't a letter.
% C_TEXT = CLEANTEXT(T,SUBSTITUTE) will replace the characters in the T
%   with characters specified in SUBSTITUTE.  SUBSTITUTE should be a cell
%   array.  The first element in each row is a string containing
%   the characters to be replaced.  The second element in each row is the
%   character used as the replacement.  If the second element is an empty
%   string, then the original characters are removed.
%
%  Example:
%
%  t = 'invalidname.123<>   ';
%  ct = cleantext(t);
%  ct2 = cleantext(t,{'<>.','_'});

if nargin <2,
    subs = {'~`!@#$%^&*()-=+|\{}[]:;''"<>,./? ','_'};
end


% remove all whitespace characters.
txt(isspace(txt)) = '';

for outer = 1:size(subs,1),
    if numel(subs{outer,2}) > 1,
        error('Each replacement should be a single character.')
    end

    for inner = 1:length(subs{outer,1}),
        if isempty(subs{outer,2})
            txt(txt == subs{outer,1}(inner)) = '';
        else
            txt(txt == subs{outer,1}(inner)) = subs{outer,2};
        end
    end
end


% if the first character isn't a letter, prepend 'var_'
%if ~isletter(txt(1))
%    txt = ['var_' txt];
%end
