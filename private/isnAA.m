function [no] = isnAA(x)
no=1;
% if (isa(x,'uint8')), 
if (isreal(x)), 
	if (x>=1&x<=23)
		no=0;
	end
end