function [t,codon] = codonise61(s,icode)
%CODONISE61 - codonises seq into 61-encoding seq (e.g., [1 1 1]->1, [4 4 4]->61)
%
%See also: DECODONISE61

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2), icode=1; end

if (ischar(s)), s=i_encode_n(s); end

[TABLE,CODON] = codontable;
stops=find(TABLE(icode,:)=='*');

t=codonise64(s);
[n,m]=size(t);
for (i=1:n),
for (j=1:m),
   %fprintf('%d %d, %d-%d\n',i,j,t(i,j), sum(t(i,j)>stops));
   t(i,j) = t(i,j)-sum(t(i,j)>stops);
end
end

if (nargout>1)
	codon={'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT'...
	      'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT'...
	      'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT'...
	      'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};
	codon(stops)=[];
end
