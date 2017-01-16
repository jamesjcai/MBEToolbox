function [S,n]=genpept2cds(peptid)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

%peptid='NP_055395';
%peptid='NP_817126';
a=getgenpept(peptid);

x=a.DBSource;
s=find(x==' ');
mrnaid=x(s(end)+1:end-2);


b=getgenbank(mrnaid);
seq=b.Sequence;
loc=b.CDS.indices;

S='';
n=size(loc,1);
for k=1:n
    S=strcat(S,seq(loc(k,1):loc(k,2)));
end

