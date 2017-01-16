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

