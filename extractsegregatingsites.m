function [aln2]=extractsegregatingsites(aln,biallelic)
%EXTRACTSEGREGATINGSITES - Extract segregating sites
%A site is considered a segregating site if there are two or more nucleotides
%at that site in a comparison of m sequences.
%
% Syntax: [aln2]=extractsegregatingSites(aln)
%
% Inputs:
%    aln     - Alignment structure
%
% Outputs:
%    aln2     - New alignment including segregating sites only
%
% See also: EXTRACTINVARIABLESITES

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if nargin<2, biallelic=0; end

if (isstruct(aln)),
%	if ~(isfield(aln,'seqtype'))
%            aln2.seqtype=1;
%	else
%	      aln2.seqtype = aln.seqtype;
%	end
%	aln2.geneticcode = 0;
%	aln2.seqnames =	 aln.seqnames;

    aln2=aln;
	[n,m] = size(aln.seq);
	aln2.seq=zeros(n,0);
	k=0;
	for j=1:m
		minnt = min(aln.seq(:,j));
		maxnt = max(aln.seq(:,j));
		if (minnt>0 & maxnt < 5)
			if (minnt~=maxnt)
			k = k+1;
			aln2.seq(:,k) = aln.seq(:,j);
            end
		end
    end

    if biallelic
        aln2.seq=i_onlybiallelic(aln2.seq);
    end
else

    seq=aln;
	[n,m] = size(seq);
	seq2=zeros(n,0);

	k = 0;
	for j=1:m
		minnt = min(seq(:,j));
		maxnt = max(seq(:,j));
		if (minnt>0 & maxnt < 5)
			if (minnt~=maxnt)
			k = k+1;
			seq2(:,k) = seq(:,j);
            end
		end
    end

    if biallelic
        [seq2]=i_onlybiallelic(seq2);
    end

    aln2=seq2;
end


function [seq2]=i_onlybiallelic(seq2)
        idx=[];
        for (k=1:size(seq2,2))
            if (length(unique(seq2(:,k)))>2)
                idx=[idx,k];
            end
        end
        seq2(:,idx)=[];






