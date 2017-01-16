function [s] = decodeseq(s0,seqtype)
%ENCODESEQ - Converts nucleotide from letters to integer.
%
% Syntax: [s] = encodeseq(s0,seqtype)
%
% Inputs:
%    s0     - Sequences letter representation
%
% Outputs:
%    s    - Sequences integer representation
%
%
% See also: CODONISESEQ, ENCODEALN

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2)
    seqtype=1;
end

if (isstruct(s0)),
      s0=s0.seq;
      seqtype=s0.seqtype;
end

[NT,AA] = seqcode;

switch (seqtype)
    case (1)
	s=NT(s0);
    case (2)
	s=NT(s0);
    case (3)
	s=AA(s0);
end




