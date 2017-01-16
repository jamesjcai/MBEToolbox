function t=codonise64(s)
%CODONISESEQ - Codonises sequence(s)
%
% Syntax:  t=codonise64(s)
%
% Inputs:
%    s   - Input sequence(s)
%
% Outputs:
%    t   - Codonised sequence(s)
%
% See also: ENCODESEQ, ENCODEALN

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

[n,m]=size(s);
if (mod(m,3)>0) error('length of coding sequence cannot divide by 3!'); end
t = zeros(n,m/3);

p=1:3:m; q=1:m/3;
t(:,q)=(s(:,p)-1).*16+(s(:,p+1)-1).*4+(s(:,p+2)-1)+1;
t(t>64)=65;