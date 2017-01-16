function plotsimilarity(aln)
%PLOTSIMILARITY - plots the running average of the similarity among the sequences in a multiple sequence alignment.
%PlotSimilarity calculates the average similarity among all members of a group 
%of aligned sequences at each position in the alignment, using a user-specified 
%sliding window of comparison. The window of comparison is moved along all 
%sequences, one position at a time, and the average similarity over the entire 
%window is plotted at the middle position of the window. 

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

slidingfun(@i_similarity,aln);


function [d] = i_similarity(aln)
[d] = pairavg(@i_dn_jc_pair,aln);
%%Another way
%if (isstruct(aln)), seq=aln.seq; else seq=aln; end
%[n,m]=size(seq);
%[D] = triu(dn_jc(seq));
%d=sum(D(:))/nchoosek(n,2);


function [res] = i_dn_jc_pair(seq)
d=dn_jc(seq);
res=d(1,2);