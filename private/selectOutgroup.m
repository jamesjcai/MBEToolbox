function [outgrno]=selectOutgroup(aln)

if (isstruct(aln)), 
	seq=aln.seq; 
	seqnames=aln.seqnames;
else 
	seq=aln;
	[n,m]=size(seq);
	for (k=1:n),
		seqnames(k)=int2str(k);
	end
end
outgrno=listdlg('ListString',seqnames,...
                'SelectionMode','single',...
                'PromptString','Select an outgroup:');
