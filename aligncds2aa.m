function [aln] = aligncds2aa(AlnNT,AlnAA)
%aligncds2aa - [aln] = aligncds2aa(AlnNT,AlnAA)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(AlnAA.seq);
aln=copyalnheader(AlnNT);
S=ones(n,m*3)*5;

SA = AlnAA.seq;
SN = AlnNT.seq;

aagap=i_getcode4gap('PROTEIN');

for (i=1:n),
x=1;
y=1;
for (j=1:m),
	if ~(SA(i,j)==aagap)
		S(i,x)=SN(i,y);
		S(i,x+1)=SN(i,y+1);
		S(i,x+2)=SN(i,y+2);
		x=x+3;
		y=y+3;
	else
		x=x+3;
	end
end
end
aln.seq = S;