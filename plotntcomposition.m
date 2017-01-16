function plotntcomposition(aln,type)
%PLOTNTCOMPOSITION - Plot nucleotide composition
%
% Syntax: plotntcomposition(aln,type)
%
% Inputs:
%    aln      - Alignment structure
%    type     - Style option
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


if (isstruct(aln)),
      seq=aln.seq;
else
      seq=aln;
end

if (nargin<2), type=1; end

switch (type)
    case (1)
	[N,Nf,NN,NNf] = ntcomposition(seq(1,:));
	bar([1 2 3 4], Nf, 0.4);
	title('Nucleotide composition');
	set(gca,'XTickLabel',{'A','C','G','T'});

	xlabel('Nucleotide (A=1, C=2, G=3, T=4)'); ylabel('Frequency (%)');
    case (2)
	info={'All positions' '1st positions' '2nd postitions' '3rd positions'};
	[N,Nf,NN,NNf] = ntcomposition(seq(1,:));
	subplot(1,4,1), bar([1 2 3 4], Nf, 0.4)
	title(info(1));
	set(gca,'XTickLabel',{'A','C','G','T'});
	ylabel('Frequency (%)');
	% xlabel(addlabel);
	for (k=1:3)
		seq2=extractpos(seq,k);
		[N,Nf,NN,NNf] = ntcomposition(seq2(1,:));
		subplot(1,4,k+1), bar([1 2 3 4], Nf, 0.4)
		set(gca,'XTickLabel',{'A','C','G','T'});
		title(info(k+1));
	end
end

