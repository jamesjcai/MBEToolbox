function plotdiplomo(aln)
%PLOTDIPLOMO - Scatter plot comparing different distance measures with each other
%DIPLOMO (DIstance PLOt MOnitor) compares different distance measures with each
%other by displaying them as a scatter plot. Originally it is an interactive
%computer graphics application, written by Georg Weiller, of the Bioinformatics
%Laboratory, Australian National University, Canberra, Australia
%(weiller@rsbs-central.anu.edu.au), uses the distance plot principle in a new
%and powerful approach to phylogenetic analysis.
%
% Syntax: plotdiplomo(aln)
%
% Inputs:
%    aln      - Alignment structure
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


[Aln1]=extractpos(aln,3);
D1 = dn_pdist(Aln1);

[aln2] = translatealn(aln);
D2 = dp_pdist(aln2);

[n,m]=size(aln.seq);
k=1;
for (i=1:n)
    for (j=i:n)
	if (i ~= j)
		prt(k,1)=D1(i,j);
		prt(k,2)=D2(i,j);
	        k=k+1;
	end
    end
end


plot(prt(:,1), prt(:,2),'+','MarkerSize',5);
info = 'Scatter plot comparing two different distances';
title(info);
xlabel('Nucleotide (3rd codon position) p-distance'); ylabel('Amino acid p-distance');