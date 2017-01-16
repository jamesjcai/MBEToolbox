function [p0,q0,p2,q2,p4,q4,m0,m2,m4] = countdegeneratesites(S)
%COUNTDEGENERATESITES - Counts degenerate sites in two aligned DNA sequences.
%This function used by DC_LI93 and DC_LI85
% Syntax: [p0,q0,p2,q2,p4,q4,m0,m2,m4] = countdegeneratesites(S)
%
% Inputs:
%    S    - Matrix of two sequences
%
% Outputs:
%    p0   - Transitions at 0-fold degenerate sites
%    q0   - Transversions at 0-fold degenerate sites
%    p2   - Transitions at 2-fold degenerate sites
%    q2   - Transversions at 2-fold degenerate sites
%    p4   - Transitions at 4-fold degenerate sites
%    q4   - Transversions at 4-fold degenerate sites
%
% See also: DC_PAMILO_BIANCHI_LI93, DC_LI_WU_LUO85

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (isstruct(S)), S=aln.seq; end

[n,m]=size(S);
if (n~=2) error('This function requires the number of sequences is 2'); end
S2=codonise64(S);
marker=zeros(n,m);

% degtable corresponding to the following codon matrix:
% codonChar = {'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA'
% 'ATC' 'ATG' 'ATT' 'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA'
% 'CTC' 'CTG' 'CTT' 'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA'
% 'GTC' 'GTG' 'GTT' 'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA'
% 'TTC' 'TTG' 'TTT'};

degtable=[0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 2 0 2; 0 0 2;
2 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 0; 0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 2;
0 0 4; 0 0 4; 0 0 4; 0 0 4; 2 0 4; 0 0 4; 2 0 4; 0 0 4; 2 0 4; 0 0 4;
2 0 4; 0 0 4; 0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4;
0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 0; 0 0 2;
0 0 0; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 0; 0 0 2; 0 0 0; 0 0 2;
2 0 2; 0 0 2; 2 0 2; 0 0 2; 0 0 0];

TS = [0, 0, 1, 0;
      0, 0, 0, 1;
      1, 0, 0, 0;
      0, 1, 0, 0];

TV = [0, 1, 0, 1;
      1, 0, 1, 0;
      0, 1, 0, 1;
      1, 0, 1, 0];

for (p=1:n),
	for (q=1:m/3)
	deg = degtable(S2(p,q),:);
	q=(q-1)*3+1;	% jump index p every next three bases.
	marker(p,q)=deg(1,1); marker(p,q+1)=deg(1,2); marker(p,q+2)=deg(1,3);
	end
end

pick0=marker==0;
pick2=marker==2;
pick4=marker==4;

m0=sum(sum(pick0));
m2=sum(sum(pick2));
m4=sum(sum(pick4));


tsmaker=zeros(1,m);
tvmaker=zeros(1,m);

for (i=1:m)
	tsmaker(i)=TS(S(1,i),S(2,i));
	tvmaker(i)=TV(S(1,i),S(2,i));
end


p0 = sum(pick0(1,:).*tsmaker) + sum(pick0(2,:).*tsmaker);
q0 = sum(pick0(1,:).*tvmaker) + sum(pick0(2,:).*tvmaker);

p2 = sum(pick2(1,:).*tsmaker) + sum(pick2(2,:).*tsmaker);
q2 = sum(pick2(1,:).*tvmaker) + sum(pick2(2,:).*tvmaker);

p4 = sum(pick4(1,:).*tsmaker) + sum(pick4(2,:).*tsmaker);
q4 = sum(pick4(1,:).*tvmaker) + sum(pick4(2,:).*tvmaker);


p0=p0/m0; q0=q0/m0; p2=p2/m2; q2=q2/m2; p4=p4/m4; q4=q4/m4;

