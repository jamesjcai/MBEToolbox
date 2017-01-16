function [TABLE,CODON] = codontable(idx)
%CODONTABLE - Return codon matrix and genetic tables matrix
%
% Syntax: [TABLE,CODON] = codontable
%
% See also: GETSEQCODE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $




 TABLE = ['KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF-';
          'KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF-';
          'NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF-';
          'KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF-';
          'NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF-';
          'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF-'];

if nargin>0
    TABLE=TABLE(idx,:);
end

if (nargout>1)
	CODON=[1 1 1; 1 1 2; 1 1 3; 1 1 4; 1 2 1; 1 2 2; 1 2 3; 1 2 4;
	       1 3 1; 1 3 2; 1 3 3; 1 3 4; 1 4 1; 1 4 2; 1 4 3; 1 4 4;
	       2 1 1; 2 1 2; 2 1 3; 2 1 4; 2 2 1; 2 2 2; 2 2 3; 2 2 4;
	       2 3 1; 2 3 2; 2 3 3; 2 3 4; 2 4 1; 2 4 2; 2 4 3; 2 4 4;
	       3 1 1; 3 1 2; 3 1 3; 3 1 4; 3 2 1; 3 2 2; 3 2 3; 3 2 4;
	       3 3 1; 3 3 2; 3 3 3; 3 3 4; 3 4 1; 3 4 2; 3 4 3; 3 4 4;
	       4 1 1; 4 1 2; 4 1 3; 4 1 4; 4 2 1; 4 2 2; 4 2 3; 4 2 4;
	       4 3 1; 4 3 2; 4 3 3; 4 3 4; 4 4 1; 4 4 2; 4 4 3; 4 4 4];
end




% GENETICTABLE = ['K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' '*' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' '*' 'S' '*' 'S' 'M' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'M' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'T' 'T' 'T' 'T' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'S' 'S' 'S' 'S' 'M' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' 'Q' 'Y' 'Q' 'Y' 'S' 'S' 'S' 'S' '*' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'N' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'S' 'S' 'S' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'C' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' '*' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'S' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' '*' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'G' 'S' 'G' 'S' 'M' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'N' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'S' 'S' 'S' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' 'Y' 'Y' '*' 'Y' 'S' 'S' 'S' 'S' 'W' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-';
% 'K' 'N' 'K' 'N' 'T' 'T' 'T' 'T' 'R' 'S' 'R' 'S' 'I' 'I' 'M' 'I' 'Q' 'H' 'Q' 'H' 'P' 'P' 'P' 'P' 'R' 'R' 'R' 'R' 'L' 'L' 'L' 'L' 'E' 'D' 'E' 'D' 'A' 'A' 'A' 'A' 'G' 'G' 'G' 'G' 'V' 'V' 'V' 'V' '*' 'Y' 'Q' 'Y' 'S' 'S' 'S' 'S' '*' 'C' 'W' 'C' 'L' 'F' 'L' 'F' '-'];
% CD = {'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT'...
%       'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT'...
%	'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT'...
%	'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};