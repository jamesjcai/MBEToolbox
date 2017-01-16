function [s] = encodeseq(s0,seqtype)
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
s=s0;
switch (seqtype)
    case (1)
         s = i_encode_n(s0);
    case (2)
         s = i_encode_n(s0);
    case (3)
         s = i_encode_a(s0);
end