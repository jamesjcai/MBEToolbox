function [D,V]=dn_pdist(aln,gappenalty)
%DN_PDIST - p-distances (nucleotide)
%
% Syntax: [D,V]=dn_pdist(aln,gappenalty)
%
% Inputs:
%    aln          - Alignment structure
%    gappenalty   - (optional) Gap penalty
%
% Outputs:
%    D      - Distance matrix
%
% See also: DP_PDIST

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[n,m] = size(seq);
if(nargin==1) gappenalty=0; end

[N,v,gap] = dn_ntdiff(aln);
D = N./((m-gap)+(gap.*gappenalty));
if (nargout>1),
	V=D.*(1-D)./m;
end