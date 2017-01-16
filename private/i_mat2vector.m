function [s] = i_mat2vector(S)

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

s=[];
[n,m]=size(S);
for (k=1:n),
      s=cat(2,s,removeblanks(S(k,:)));
end