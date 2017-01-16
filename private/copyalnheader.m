function [aln2] = copyalnheader(aln)
%COPYALNHEADER - Copy header information of an alignment structure into another.
%
% Syntax: [aln2] = copyalnheader(aln)
%
% Inputs:
%    aln      - Alignment structure
%
% Outputs:
%    aln2     - Alignment structure
%
%
% See also: 

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/28/2005

aln2.seqtype=aln.seqtype;
aln2.geneticcode = aln.geneticcode;
aln2.seqnames = aln.seqnames;