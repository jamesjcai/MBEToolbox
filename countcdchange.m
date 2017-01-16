function [D,gap]=countcdchange(s1,s2)
%COUNTCDCHANGE - Counts codon (61 sense codons) replacements between two sequences
%
% Syntax: [D,gap]=countcdchange(s1,s2)
%
% Inputs:
%    s1   - Sequence 1 vector
%    s2   - Sequence 2 vector
%
% Outputs:
%    D    - Matrix of codon difference
%    gap  - Number of gaps
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
	[D,gap]=countchange(s1,s2,61);
else
	[D]=countchange(s1,s2,61);
end
