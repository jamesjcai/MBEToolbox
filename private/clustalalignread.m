function S = clustalalignread(filename)
% CLUSTALALIGNREAD - imports files containing Clustal alignments
% S = CLUSTALALIGNREAD(ALIGNFILE) will import the file ALIGNFILE and
%   convert it to a structure S.  ALIGNFILE should be in the Clustal
%   alignment format.  TS will contain one field for each sequence 
%   in the file, and one for the consensus line.

txt = textread(filename,'%s','delimiter','\n','whitespace','');

% check header line, confirm that it's a ClustalW file
if ~strncmpi(txt{1},'CLUSTAL',7)
    error('Header does not match CLUSTAL format')
end

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

% the first consensus line is in mt(1)-1
cons_loc = mt(1)-1;

% there are cons_loc-1 sequences
num_seq = cons_loc-1;

% create empty structure
S = struct;
for s = 1:num_seq,
    % make the name into a MATLAB-acceptable variable name
    name = cleantext(strtok(txt{s}),{'|','_';'.',''});
    
    % initialize field to hold sequence
    S.(name) = '';
    for r = s:cons_loc:size(txt,1),
        % make sure that there aren't sequence numbers at the end
        S.(name) = [S.(name) deblank(strtok(txt{r}(17:end),'0123456789'))];      
    end
end

% consensus line
S.consensus = '';
for r = cons_loc:cons_loc:size(txt,1)
    S.consensus = [S.consensus txt{r}(33:end)];
end






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
if ~isletter(txt(1))
    txt = ['var_' txt];
end








