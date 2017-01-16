function [aln2,sites]=extractinformativesites(aln,variants)
%EXTRACTINFORMATIVESITES - Extract informative sites
%Sites returns the indices of columns with at least 2 bases appearing at least
%twice each. These are the informative sites of the method of maximum parsimony.
%
% Syntax: [aln2]=extractinformativesites(aln)
%
% Inputs:
%    aln     - Alignment structure
%
% Outputs:
%    aln2     - New alignment including informative sites only
%
% See also: EXTRACTSEGREGATINGSITES

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2), variants=2; end
if (variants<2), variants=2; end

%if ~(isvalidaln(aln,'CODING')|isvalidaln(aln,'NONCODING'))
%	error ('ERROR: Not nucleotide seq')
%end
% if (hasgap(aln))
% 	error ('ERROR: Sequences cannot contain gaps. Gaps can be removed by using REMOVEGAPS.')
% end

seqarray=aln.seq;
[m,n]=size(seqarray);
if m<4
   disp('Too few sequences.')
end

aln2=aln;
%aln2.seqtype = aln.seqtype;
%aln2.geneticcode = 0;
%aln2.seqnames =	aln.seqnames;

% informative.m
%
% usage: sites=informative(seqarray)
%
% Each row of seqarray should be a sequence in A,G,C, and T,
% with one row per taxa;
% sites returns the indices of columns with at least 2 bases
% appearing at least twice each. These are the informative
% sites of the method of maximum parsimony.
%
% 8/2/03

%if variants=2, two variants;
%if variants=3, three variants;
%if variants=4, four variants;

variants=variants-1;

As=( (sum(seqarray==1)) > variants );	% A
Gs=( (sum(seqarray==3)) > variants );	% G
Cs=( (sum(seqarray==2)) > variants );	% C
Ts=( (sum(seqarray==4)) > variants );	% T

sites=find( (As+Gs+Cs+Ts) > 1 );
aln2.seq=aln.seq(:,sites);