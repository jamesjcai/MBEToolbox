function [t]=codonise64(s)
%CODONISE64 - Codonises sequence(s)
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

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[n,m]=size(s);
if (mod(m,3)>0), error('length of coding sequence cannot divide by 3!'); end
t = zeros(n,m/3);

p=1:3:m; q=1:m/3;
t(:,q)=(s(:,p)-1).*16+(s(:,p+1)-1).*4+(s(:,p+2)-1)+1;
%t(unique(ceil(find(s>4)./3)))=65;
t(t>64)=65;
