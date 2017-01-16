function [D,V]=dp_poisson(aln)
%DP_POISSON - Poisson Correction (PC) distance (amino acid)
% Poisson Correction (PC) distance
% The Poisson correction distance corrects for multiple substitutions
% at the same site, but it assumes equality of substitution rates
% among sites and assumes equal amino acid frequencies.
%
% Syntax: [D,V]=dp_poisson(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    D     - Distance matrix
%    V     -
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
% DVAR = zeros(n);


[P,v]=dp_pdist(aln);
D= -log(1-P);
V=P./((1-P)*m);

% SE=sqrt(V);