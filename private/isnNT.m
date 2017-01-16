function [no] = isnNT(n)
no=1;
% if (isa(n,'uint8')), 
if (isreal(n)), 
	if (n>=1&n<=5)
		no=0;
	end
end