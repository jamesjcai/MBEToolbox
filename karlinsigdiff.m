function [d] = karlinsigdiff(s1,s2)
%KARLINSIGDIFF - Return the difference of Karlin genomic signatures between two given sequences
%Difference of Karlin genomic signatures is a measure of the genomic difference
%between two sequences fand g (from different organisms or from different
%regions of the same genome). It is calculated as the dinucleotide average
%absolute relative abundance difference:
%  delta*(f,g) = (1/16)sum(|rho*(f)-rho*(g)|)
%where the sum extends over all dinucleotides.
%
%For details see "Global dinucleotide signatures and analysis of genomic
%heterogeneity" Samuel Karlin - Current Opinion in Microbiology 1998,
%1:598-610.
%
% Syntax: [d] = karlinsigdiff(s1,s2)
%
% Inputs:
%    s1    - Nucleotide sequence 1
%    s2    - Nucleotide sequence 2
%
% Outputs:
%    d    - Dinucleotide average absolute relative abundance difference
%
% See also: CALCULATEKARLINSIG

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


K1 = karlinsig(s1);
K2 = karlinsig(s2);
d=sum(sum(abs(K1-K2)))/16;