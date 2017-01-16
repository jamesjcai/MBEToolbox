%% MBEToolbox DEMO - Basic sequence statistics
% Welcome to MBEToolbox.  This is a demonstration of
% MBEToolbox functions for basic sequence manipulation and
% statistics.
%
% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



cd(fileparts(which('mbe_demo2')))
% Current working directory will be changed.

%%
% First, let's read a FASTA formatted file into a MATLAB structure
% represents an alignment and view the sequences in this alignment.

filename=fullfile('seq_examples','HCV_aligned.fas');
aln = readfasta(filename,2,1)
viewseq(aln)

%%
% Now let's view some operations working on single sequence.
%
% Notice we use the following command to extract the first
% sequence in alignment.

seq = aln.seq(1,:)

%%
% We count its nucleotide and dinucleotide composition.

[N,Nf,NN,NNf] = ntcomposition(seq)

%%
% Now let's plot nucleotide composition in a bar chart
%
% Notice this plot working on the first sequence in alignment.

plotntcomposition(aln,2)

%%
% Karlin genomic signatures for a given sequence
%
% K - A 4x4 matrix dinucleotide relative abundance values.

[K] = karlinsig(seq)

%%
% Performs sliding window analysis of GC deviation

subplot(1,1,1);
plotslidingwin(aln,'GCDeviation',50)

%%
% Codon volatility can be calculated on single coding sequence or
% a set of sequences all together

[v,V] = codonvolatility(aln)


%%
% Return alignment with rev. com. strand of nucleotide sequences
aln2 = revcomseq(aln)
viewseq(aln2)

%%
% Extract coding position 3

aln2 = extractpos(aln,3)
viewseq(aln2)


%%
% Counts GC content at fourfold degenerate sites

[V,P] = gc4(aln)

%%
% Translate coding DNA sequence into protein sequence in an alignment

aln2 = translatealn(aln)
viewseq(aln2)

%%
% Thank you for viewing this introduction to MBEToolbox basic sequence
% manipulation and statistics.
