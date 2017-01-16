function [aln] = arrseq2tre(aln,tr)
%Arrange seqences in the order of nodes on tree
% [aln] = arrseq2tre(aln,tree)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


[aln.seqnames,idx]=sort(aln.seqnames);
aln.seq=aln.seq(idx,:);

names={};
if (nargin>1)
	[dummy,dummy,dummy,nm]=parsetree(tr);
	if ~(iscell(nm)),
	    for (k=1:size(nm,1)),
	          names(k)={deblank(nm(k,:))};
	    end
	end
	[dummy,idx]=sort(names);

	

	% now rearrange seq and seqnames according to idx, index of tree nodes.
	for (k=1:length(idx)),
	if ~(ismember(names(k),aln.seqnames)),
	      error('Not match')
	end

	% now rearrange seq and seqnames according to idx, index of tree nodes.
		for (k=1:length(idx)),
			xnames(idx(k))={char(aln.seqnames(k))};
			xseq(idx(k),:)=aln.seq(k,:);
		end
	end
	

	aln.seq=xseq;
	aln.seqnames=xnames;
end