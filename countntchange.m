function [D,gap]=countntchange(s1,s2)
%COUNTNTCHANGE - Count nucleotide changes in two DNA sequences
%D is a 4x4 array, with bases in seq1 along top, seq2 along side,
%in order A,C,G,T.
%
% Syntax: [D,gap]=countntchange(s1,s2)
%
% Inputs:
%    s1      - Sequence 1 vector
%    s2      - Sequence 2 vector
%
% Outputs:
%    D       - Codon Adaptation Index value
%    gap     - Codon Adaptation Index value
%
%
% See also: COUNTAACHANGE COUNTCDCHANGE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargout>1)
	[D,gap]=countchange(s1,s2,4);
else
	[D]=countchange(s1,s2,4);
end
