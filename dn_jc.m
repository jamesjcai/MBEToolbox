function [D,VarD]=dn_jc(aln,gappenalty,alpha)
%DN_JC - Jukes-Cantor Distance
%The Jukes-Cantor correction is used in nucleotide distances and the synonymous and nonsynonymous
%substitution distances.  If the proportion of sites that are different (nucleotides, synonymous,
%or nonsynonymous) is greater than or equal to 75%, the Jukes-Cantor correction cannot be applied.
%This seems to have happened for some pairs in your data.  If you wish to know which pair(s) of
%sequences has this problem, please use the Distances|Pairwise option.  All such pairs will be
%shown in the Distance Matrix with a NaN.
%
% Syntax: [D,VarD]=dn_jc(aln,gappenalty,alpha)
%
% Inputs:
%    aln          - Alignment structure
%    gappenalty   - (optional) Gap penalty
%    alpha        - (optional) Shape parameter of gamma distribution
%
% Outputs:
%    D      - Distance matrix
%    VarD   - Variance of distance
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


if(nargin==1)
	[P,v]=dn_pdist(aln);
	W1=1-4*P./3;
	D=(-3/4)*i_safelog(W1);
elseif(nargin==2)
	[P,v]=dn_pdist(aln,gappenalty);
	W1=1-4*P./3;
	D=(-3/4)*i_safelog(W1);
elseif(nargin==3)
	[P,v]=dn_pdist(aln,gappenalty);
	W1=1-4*P./3;
	D=(3/4)*alpha*((max(realmin,W1).^(-1/alpha)-1));
end

for k=1:size(D,1)
    D(k,k)=0;
end

if (nargout>1),
	[n,m] = size(aln.seq);
	VarD=(P.*(1-P))./(m*(W1.^2));
end
