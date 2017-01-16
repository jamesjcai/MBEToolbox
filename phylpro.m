%{
function [] = phylpro()
%PHYLPRO - Stepwise detection of recombination breakpoints using phylogenetic profiling (Phylpro) at each step
Arguments
input_file 	character string indicating the name of a Phylip format data input file
breaks 	an integer vector of ordered site(s) just before the previously declared breakpoints
winHalfWidth 	the window half width to use
permReps 	the number of Monte Carlo replicates to use for the permutation distribution

The phylpro function implements phylogenetic profiling (Phylpro) for detecting recombination breakpoints (Weiller 1998) using a moving window of fixed width. Breakpoints detected in previous steps of a stepwise search may be conditioned upon.

For a given position of the moving window on the sequence alignment, and for a given “target” sequence, a correlation is computed to compare two distance vectors: the distance between the target sequence and all other sequences in the left half-window and the distance between the target sequence and all others in the right half-window. The pair-wise distance measure used is the proportion of sites at which the sequences differ. Discordance between the two distance vectors may reflect a recombination event, located at the window centre, in the history of the target sequence. The minimum correlation over all target sequences is regarded as a summary of the evidence for recombination at the window centre. The individual correlations for the target sequences may also be of interest for suggesting sequence segments that descend from historical recombination events. Significance of observed correlation statistics is assessed by a Monte Carlo permutation test. When conditioning on breakpoints proposed at previous steps of a stepwise search, permutation is restricted to sites within blocks defined by the previously proposed breakpoints, as described by Graham et al. (2004).
Value
polyposn 	The site numbers of all ungapped polymorphic sites in the alignment
corrs 	Observed correlations that exceed the 90th percentile of the permutation null distribution
winlocs 	Window centres corresponding to the correlations in corrs
target.seqs 	The target sequence that lead to a significant correlation in corrs
quants 	90th and 95th percentiles of the permutation distribution
Author(s)

Brad McNeney <mcneney@stat.sfu.ca>, Jinko Graham <jgraham@stat.sfu.ca>, Sigal Blay <sblay@sfu.ca>
References

Graham J, McNeney B and Seillier-Moiseiwitsch F (2004). Stepwise detection of recombination breakpoints in sequence alignments. Bioinformatics Sep 23; [Epub ahead of print]

Weiller G (1998). Phylogenetic profiles: A graphical method for detecting genetic recombination in homologous sequences. Mol Biol Evol, 15:326-335.

http://stat-db.stat.sfu.ca/stepwise 
%}