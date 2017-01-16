function [aln2] = translatealn(aln)
%TRANSLATEALN - Translate coding DNA sequence into protein sequence in an alignment
%
% Syntax: [aln2] = translatealn(aln)
%
% Inputs:
%    aln      - Alignment structure (protein-coding nucleotide sequences)
%
% Outputs:
%    aln2     - New alignment structure (protein sequences)
%
% See also: REVERSEALN, REVCOMALN

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



%if ~(isvalidaln(aln,'CODING'))
%	error ('ERROR: Not coding seq')
%end

aln2=aln;
aln2.seqtype=3;
aln2.seq = translateseq(aln.seq,aln.geneticcode);
