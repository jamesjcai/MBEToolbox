function [D,VarD]=dn_tajima_nei84(aln)
%DN_TAJIMA_NEI84 - Tajima & Nei 84 Distance
%In the real data, nucleotide frequencies often deviate substantially from
%0.25. In this case the Tajima-Nei distance (Tajima and Nei 1984) gives a
%better estimate of the number of nucleotide substitutions than the Jukes-
%Cantor distance. Note that it assumes equality of substitution rates among
%sites and between transitionalal and transversional substitutions. This
%method is only for nucleotide sequences.
%
% Syntax: [D,VarD]=dn_tajima_nei84(aln)
%
% Inputs:
%    aln    - Alignment structure
%
% Outputs:
%    D      - Distance matrix
%    VarD   - Variance of distance
%
% REF: Li WH (1997) Molecular Evoluton p. 81
%      Nei and Kumar (2000), page 38.
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
if (n<2), error('At least two sequences.'); end

D=zeros(n,n);
VarD=zeros(n,n);


for i=1:n-1
for j=i+1:n
	S1=S(i,:); S2=S(j,:);
	[d,v] = i_SeqPairDistanceTajimaNei84(S1,S2);
	D(i,j)=d; D(j,i)=d;
	VarD(i,j)=v; VarD(j,i)=v;
end
end







function [d,v] = i_SeqPairDistanceTajimaNei84(S1,S2)

[n,m] = size(S1);

q=zeros(1,4);
% -------- y(i,j), y(j,i)
Y=zeros(4,4);
for (i=1:4),
	q(i) = (sum(S1==i)+sum(S2==i))/(length(S1)+length(S2));	% Estimated equilibrium frequence of A,C,G,T
	for (j=1:4),
		Y(i,j) = length(find(S1==i&S2==j));
	end
end
% -------- x(i,j)
X = triu(Y)+tril(Y)';

p=sum(S1~=S2)/m;

% -------- h
h=0;
for (i=1:4),
for (j=i:4),
	h=h+(X(i,j)*X(i,j))/(2*q(i)*q(j));
end
end

% -------- b1
b1=1-sum(q.^2);
b=(b1+(p^2)/h)/2;

d=-1*b*i_safelog(1-p/b);
v=b^2*p*(1-p)/((b-p)^2*m);
