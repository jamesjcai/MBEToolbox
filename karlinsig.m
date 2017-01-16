function [Kf,Kn] = karlinsig(S)
%KARLINSIG - Returns the Karlin genomic signatures for a given sequence
%Karlin's genomic signatures measure dinucleotide bias (dinucleotide relative
%abundance values) of a DNA sequence. They are assessed through the odds ratos
%  rho(XY) = f(XY)/f(X)f(Y)
%where f(XY) is the frequency of the dinucleotide XY and f(X) is the frequency
%of the nucleotide X. For double-stranded DNA sequences, a symmetrized version
%{rho*(XY)} is computed from corresponding frequencies of the sequence
%concatenated with its inverted complementary sequence.
%
%For details see "Global dinucleotide signatures and analysis of genomic
%heterogeneity" Samuel Karlin - Current Opinion in Microbiology 1998,
%1:598-610.
%
%This function DOES NOT contain the code for calculating the absolute "genomic
%signature" difference.
%
% Syntax: [K] = karlinsig(S)
%
% Inputs:
%    S    - Nucleotide sequence
%
% Outputs:
%    K    - A 4x4 matrix dinucleotide relative abundance values.
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



% If the start/end pair doesn't give a multiple of three, bases end is moved down
% so that it is a multiple of three.
% A=1; B = B-mod((B-A+1),3);

ws = warning('off','MATLAB:divideByZero');

K=zeros(4);

S(find(S>4|S<1))=[];    %remove gaps
S=cat(2,S,revcomseq(S));
[n,m] = size(S);
m=m-mod(m,2);

for (x=1:2:m),
      K(S(1,x),S(1,x+1))=K(S(1,x),S(1,x+1))+1;
end
Kn=K;
K=K./(m/2);


F=[sum(S==1),sum(S==2),sum(S==3),sum(S==4)];
F=F./sum(F);


FF=zeros(4);
for (i=1:4),
for (j=1:4),
	FF(i,j)=F(i)*F(j);
end
end

Kf=K./FF;

warning(ws);

