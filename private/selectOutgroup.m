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


items(1).name = 'OUTGROUP:  ';
items(1).default = 1;
items(1).indent = 0;
items(1).values = seqnames;
items(1).help = 'Please select an outgroup';

title = 'Selecte Outgroup';
msg = sprintf(['Please select sequence type and genetic code']);
out = CSEFlagDialog(items, title,msg);
if ~(isempty(out)),
	outgrno=out(1).answer;
else
	outgrno=0;
end   
