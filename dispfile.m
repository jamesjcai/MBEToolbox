function dispfile(filename)
%DISPFILE -

if ~ischar(filename)
    error('mbetoolbox:InvalidInput','Input must be a character array')
end

if ~(exist(filename,'file') || exist(fullfile(cd,filename),'file')),
%  is a valid filename ?
    error('mbetoolbox:InvalidInput','Input must be a valid file')
end



file = fopen(filename, 'r');
disp(['Reading ',filename]);
txt = textread(filename,'%s','delimiter','\n','whitespace','');
for (k=1:length(txt)),
      fprintf('%s\n',char(txt(k)));
end
