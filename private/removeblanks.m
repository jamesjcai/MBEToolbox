function [out] = removeblanks(in)
%REMOVEBLANKS - Removes space within sequence,as well as, both leading and 
%trailing blanks
[n,m] = size(in);
codes = zeros(1,m);
for i = 1:m
	for j = 1:n
		if (in(j,i) == ' ')
			codes(i) = 1;
			break;
		end
	end
end
for i = 1:m
	j = m - i + 1;
   if codes(j) == 1
	in(:,j) = [];
   end
end

[r,c] = find( (in~=0) & ~isspace(in) );
if isempty(c),
	out = in([]);
else
	out = in(:,min(c):max(c));
end