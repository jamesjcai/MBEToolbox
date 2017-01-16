function [aln2] = encodealn(aln)
%ENCODEALN - Convert nucleotide in alignment to integer.
%
% Syntax: [aln2] = encodealn(aln)
%
% Inputs:
%    aln     - Alignment structure letter representation
%
% Outputs:
%    aln2    - Alignment structure integer representation
%
%
% See also: CODONISESEQ, ENCODESEQ

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if ~(aln.seqtype),
    error('Do not know the type of sequence!');
end
aln2=aln;

switch (aln2.seqtype)
    case (1)
         aln2.seq = i_encode_n(aln2.seq);
    case (2)
         aln2.seq = i_encode_n(aln2.seq);
    case (3)
         aln2.seq = i_encode_a(aln2.seq);
end
