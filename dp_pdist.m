function [D,V]=dp_pdist(aln,gappenalty)
%DP_PDIST - p-distances (amino acid)
%
% Syntax: [D,V]=dp_pdist(aln,gappenalty)
%
% Inputs:
%    aln          - Alignment structure
%    gappenalty   - (optional) Gap penalty
%
% Outputs:
%    D      - Distance matrix
%
% See also: DN_PDIST

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



[n,m] = size(aln.seq);
if(nargin==1) gappenalty=0; end

[N,v,GAP] = dp_aadiff(aln);
D = N./((m-GAP)+(GAP.*gappenalty));
V=D.*(1-D)./m;