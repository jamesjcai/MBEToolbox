function [CodonNum,CodonFreq] = codonusage(aln,printoption)
%CODONUSAGE - Codon usage.
%D is a 4x4 array, with bases in seq1 along top, seq2 along side, in order A,C,
%G,T
%
% Syntax: [CodonNum,CodonFreq] = codonusage(aln,printoption)
%
% Inputs:
%    aln         - Alignment structure
%    printoption - Print option. 0 - not print; 1 - print
%
% Outputs:
%    CodonNum    - Number of codon
%    CodonFreq   - Frequency of codon
%
% Examples:
%
% >> [CodonNum,CodonFreq] = codonusage(aln);
% >> [CodonNum,CodonFreq] = codonusage(aln,0);
% >> [CodonNum,CodonFreq] = codonusage(aln,1);
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $

if (isstruct(aln))
    if ~(isvalidaln(aln,'CODING')) error ('Not coding seq'); end
    if(nargin==1) printoption=0; end
        S=codonise64(aln.seq);
    	Seq=aln.seq;
else
	Seq=aln;
	S=codonise64(Seq);
	printoption=0;
end


[n,m]=size(Seq);
CodonNum=zeros(n,64);
for (k=1:64)
 	CodonNum(:,k)=sum((S==k),2);
end

if(nargout>1)
    CodonFreq=CodonNum./m;
end


if (printoption)
CD = {'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG'...
      'AGT' 'ATA' 'ATC' 'ATG' 'ATT' 'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC'...
      'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT' 'GAA'...
      'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT'...
      'GTA' 'GTC' 'GTG' 'GTT' 'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG'...
      'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};

	for (k=1:n)
		fprintf(['%d = %s\n'], k, char(aln.seqnames(k)));
	end
	fprintf(['\n------------------------------------------------------\n']);
	for (i=1:64)
		fprintf(['%s\t\t'], char(CD(i)));
		for (j=1:n)
			fprintf(['%d(%#2.2f)\t\t'],CodonNum(j,i),CodonFreq(j,i)*100);
		end
		fprintf(['\n']);
	end
end
