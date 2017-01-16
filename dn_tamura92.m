function [D]=dn_tamura92(aln,alpha)
%DN_TAMURA92 - Tamura 92 (T3P) Distance
%This method is only for nucleotide sequences. This method uses transition and
%transversion rates and takes into account the deviation of GC content from the
%expected value of 50 %. Gap and ambiguous positions are ignored.
%
% P = transitions/npos
% Q = transversions/npos
%
% npos - number of positions scored
%
% GC1 = GC fraction in sequence 1
% GC2 = GC fraction in sequence 2
% C = GC1 + GC2 - 2*GC1*GC2
%
% distance = -C ln(1-P/C-Q) - 0.5(1-C) ln(1-2Q)
%
% Reference:
% K. Tamura, Mol. Biol. Evol. 1992, 9, 678.
% Tamura 3-parameter distance
%
% Syntax: [D]=dn_tamura92(aln,alpha)
%
% Inputs:
%    aln          - Alignment structure
%    alpha        - (optional) Shape parameter of gamma distribution
%
% Outputs:
%    D      - Distance matrix
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
[P,Q]=countseqpq(S);
P=P./m; Q=Q./m;

% P = transitions/npos
% Q = transversions/npos
% npos - number of positions scored
% GC1 = GC fraction in sequence 1
% GC2 = GC fraction in sequence 2
% C = GC1 + GC2 - 2*GC1*GC2
% EMBOSS

WG = ones(n,n);
GC = zeros(n,1);

countGC=sum((S==2|S==3),2); % counting GC
GC = countGC./m;


for i=1:n-1
for j=i+1:n
	WG(i,j)=GC(i,1)+GC(j,1)-2*GC(i,1)*GC(j,1);
	WG(j,i)=WG(i,j);
end
end

W1=1-P./WG-Q;
W2=1-2*Q;
% CC=1-C;

if(nargin==1)
	D=-1*WG.*i_safelog(W1)-(1/2)*(1-WG).*i_safelog(W2);
elseif(nargin==2)
	error('Not implemented yet');
	% D=(alpha/2)*((WG.*W1).^(-1/alpha)+(1-WG).*W2.^(-1/alpha)-WG-1);
end

% D = (-C.)*log(1-P./C.-Q.)-(1/2)*(1-C)*log(1-2*Q);
% D=C;
% D = (-C)*log(1-P/C-Q)-(1/2)*(1-C)*log(1-2*Q);

%gc1 = count_gc(n,m,Seq(1,:))/m;
%gc2 = count_gc(n,m,Seq(2,:))/m;
%c = gc1 + gc2 - 2*gc1*gc2;
%d = (-c)*log(1-p/c-q)-(1/2)*(1-c)*log(1-2*q);

% MEGA2
% gc = (count_gc(n,m,Seq(1,:)) + count_gc(n,m,Seq(2,:)))/(2*m);
% w0 = 2*gc*(1-gc);
% w1 = 1-p/w0-q;
% w2 = 1-2*q;
% d = (-c)*log(w1)-(1/2)*(1-w0)*log(w2);