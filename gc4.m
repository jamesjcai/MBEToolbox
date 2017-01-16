function [v,p]=gc4(aln)
%GC4 - Counts GC content at fourfold degenerate sites
%
% Syntax: [v,p]=gc4(aln)
%
% Inputs:
%    aln  - Alignment structure
%
% Outputs:
%    V    - Number of GC
%    P    - Frequency of GC
%
% See also: COUNTBASECHANGE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if isstruct(aln)
    S=aln.seq;
else
    S=aln;
end

[seq0,seq2,seq4,m0,m2,m4] = extractdegeneratesites(S);
v=sum((seq4==2|seq4==3),2);
p=v./m4;