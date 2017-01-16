function twosliding(aln,id1,id2,winsize)

if (nargin<4), winsize=20; end
name1=aln.seqnames(id1);
name2=aln.seqnames(id2);
seq1=aln.seq(id1,:);
seq2=aln.seq(id2,:);

aln2=aln;
aln2.seqnames={char(name1), char(name2)};
aln2.seq=[seq1; seq2];

subplot(2,1,1); plotslidingwinkaks(aln2,winsize,1);
title([char(name1), ' -- ' ,char(name2), ' (window size = ', num2str(winsize), ')'])
subplot(2,1,2); plotslidingwinkaks(aln2,winsize,3);
