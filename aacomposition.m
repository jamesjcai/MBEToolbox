function [n,freq] = aacomposition(s)
%AACOMPOSITION - Counts AA composition
%
% Syntax: [n,freq] = aacomposition(s)
%
% Inputs:
%    s   - Protein sequence(s)
%
% Outputs:
%    n      - 1x20 vector containg number of AAs
%    freq   - Frequencies of n
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[n,m]=size(s);
%if ~(n==1)
%	error('Must be a single sequence.')
%end

n=zeros(1,20);

for (k=1:20),
      n(1,k)=sum(sum(s==k));
end
if (nargout>1), freq=n./sum(n); end
