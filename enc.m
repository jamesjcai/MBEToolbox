function [nc]=enc(seq)

%The Effective Number of codons measure is a way of analysing how biased a gene is in terms of its codon usage. Nc values range from 61 for a gene that tends to use all codons with equal frequency to 20 for a gene that is effectively using only a single codon for each amino acid.
%There is a relationship between Nc and the base composition of a gene with genes that have more biased base compositions being expected to have lower Nc values. What is usually of interest are Nc values that are lower than might be dictated by the base composition of the gene. This might be taken as evidence that there is some kind of selective pressure on the gene to use a smaller subset of codons. This selective pressure could be translational selection for 'optimal' codons. Optimal codons are those that correspond to the major abundance tRNA for that amino acid. In such circumstances, there could be a selective pressure to use a particular codon that corresponds to this tRNA. 
%The effective number of codons (NC) (Wright 1990).
%This index is a simple measure of overall codon bias and is analogous to the effective number of alleles measure used in population genetics. Knowledge of the optimal codons or a reference set of highly expressed genes is unnecessary. Initially the homozygosity for each amino acid is estimated from the squared codon frequencies (see Wright 1990 ).
%If amino acids are rare or missing, adjustments must be made. When there are no amino acids in a synonymous family, Nc is not calculated as the gene is either too short or has extremely skewed amino acid usage (Wright 1990). An exception to this is made for genetic codes where isoleucine is the only 3-fold synonymous amino acid, and is not used in the protein gene. The reported value of Nc is always between 20 (when only one codon is effectively used for each amino acid) and 61 (when codons are used randomly). If the calculated Nc is greater than 61 (because codon usage is more evenly distributed than expected), it is adjusted to 61. 

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (isstruct(seq))
    seq=seq.seq;   % when ALN as input
end

n=size(seq,1);          % number of sequences
nc=zeros(n,1);

[CodonNum,CodonFreq] = codonusage(seq);


% n is the number of codons used in the gene for the amino acid 
% pi is the frequency (within the amino acid group) of the ith codon (ni/n)



(CodonNum.*sum(CodonFreq.^2)-1)./(CodonNum-1)



deglist=zeros(1,64);
ct=codontable(1);
for k=1:64
    if ~strcmp(ct(k),'*')
      [deglist(k)]=length(find(ct(k)==ct));
    end
end   
x=[grpstats(deglist,deglist,'mode'),grpstats(deglist,deglist,'length')];
y=x(:,2)./x(:,1);

%Designated SF types 1,2,3,4, and 6 according to their repective number of
%synonymous codons

seq1=seq(1,:);
seq1c=codonise64(seq1);

deglist(seq1c)


CodonNum(find(ct(23)==ct))
