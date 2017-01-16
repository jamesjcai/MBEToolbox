function [y]=hasgap(aln)
%HASGAP - Check if alignment contains gap
%
% Syntax: [y]=hasgap(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    y     - 1 or 0
%
%
% See also: REMOVEGAPS

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


% NT = ['A' 'C' 'G' 'T' '-'];
% AA = ['A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y' '*' '-'];
if (isstruct(aln)),
	seq=aln.seq;
	if (aln.seqtype==1|aln.seqtype==2)	% DNA/RNA
		y = max(max(seq))>4;
	elseif (aln.seqtype==3)	% PROTEIN
		y = max(max(seq))>20;
	else
	    error ('No seqtype.');
	end
else
	seq=aln;
	y = max(seq(:))>4;
end
