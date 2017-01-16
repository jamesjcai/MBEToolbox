function [p,n]=gc3(aln)
%GC3 - Counts GC content at thrid positions
%
% Syntax: [p,n]=gc3(seq)
%
% Inputs:
%    aln  - Alignment structure or sequence
%
% Outputs:
%    n    - Number of GC
%    p    - Frequency of GC
%
% See also: COUNTBASECHANGE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if isstruct(aln)
    seq=aln.seq;
else
    seq=aln;
end

seq3=extractpos(seq,3);
[Ntt,Nf,NN,NNf] = ntcomposition(seq3);
n=Ntt(2)+Ntt(3);
p=n/(sum(Ntt));
