function [D]=dn_ntfreqdist(aln)
%DN_NTFREQDIST - Euclidean distances between nucleotide frequencies
%REF: (Orti G & Meyer A 1996 MBE)
%
% Syntax: [D]=dn_ntfreqdist(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    D     - Distance matrix
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (isstruct(aln)),
	S=aln.seq;
else
	S=aln;
end

[n,m] = size(S);
D = zeros(n);


for i=1:n-1
for j=i+1:n
	X1 = i_countNtComposition(S(i,:));
	X2 = i_countNtComposition(S(j,:));
	D(i,j) = sqrt(sum((X1-X2).^2));
	D(j,i) = D(i,j);
end
end


function  [X] = i_countNtComposition(Seq)

As=Seq==1;
Cs=Seq==2;
Gs=Seq==3;
Ts=Seq==4;
NA=sum(As');
NC=sum(Cs');
NG=sum(Gs');
NT=sum(Ts');

X=cat(2,NA,NC,NG,NT);