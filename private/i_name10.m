function [name2] = i_name10(name)

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

   len=length(name);
   if (len<10)
	name(len+1:10)=' ';
   elseif (len>10)
	name = char(name(1:10));	
   end
name2 = name;