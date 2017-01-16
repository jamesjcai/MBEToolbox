function [freq] = estimatefreq(S)
%ESTIMATEFREQ - Estimates base frequencies of the given sequence(s)
%Get empirical base frequencies from the data
%
% Syntax: [freq] = estimatefreq(S)
%
% Inputs:
%    S        - Sequence pair
%
% Outputs:
%    freq - Frequencies of base, e.g. [.1 .2 .3 .4]
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



N=[sum(sum(S==1)),sum(sum(S==2)),sum(sum(S==3)),sum(sum(S==4))];
freq=N./sum(N);		% Pi(A,C,G,T);
