function comphist(s,seqtype)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (nargin<2)
      seqtype=1;
end

switch (seqtype)
    case (1)
	[N,Nf,NN,NNf] = ntcomposition(s);
	%GC=(N(2)+N(3))/(sum(N));
	%fprintf('GC - %1.2f\n', GC);
	%fprintf('A - %5d (%1.2f)\n', N(1),Nf(1));
	%fprintf('C - %5d (%1.2f)\n', N(2),Nf(2));
	%fprintf('G - %5d (%1.2f)\n', N(3),Nf(3));  
	%fprintf('T - %5d (%1.2f)\n\n', N(4),Nf(4));
	plotntcomposition(s,2);
    case (2)
	[CodonNum,CodonFreq] = codonusage(s);
	CD = {'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT' 'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT' 'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT' 'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};
	barh(CodonNum)
	set(gca,'YTick',1:64)
	set(gca,'YTickLabel',char(CD));
	set(gca,'FontSize',7);
	set(gca,'YLim',[0 64+1])
	set(gca,'YDir','reverse')
    case (3)
         [n,freq] = aacomposition(s);
PROT={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
barh(n)
	set(gca,'YTick',1:20)
	set(gca,'YTickLabel',char(PROT));
	set(gca,'FontSize',7);
	set(gca,'YLim',[0 20+1])
	set(gca,'YDir','reverse')

end
