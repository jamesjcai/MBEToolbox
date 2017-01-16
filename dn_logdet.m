function [D]=dn_logdet(aln)
%DN_LOGDET - Log-det (paralinear) distance
%The LogDet model computes the distance from the determinant of the matrix of
%co-occurrence of nucleotides in the two species, according to the formula
%
%   D  = - 1/4(loge(|F|) - 1/2loge(fA1 fC1 fG1 fT1 fA2 fC2 fG2 fT2))
%
%Where F is a matrix whose (i,j) element is the fraction of sites at which base
%i occurs in one species and base j occurs in the other. fji is the fraction of
%sites at which species i has base j. The LogDet distance cannot cope with
%ambiguity codes. It must have completely defined sequences. One limitation of
%the LogDet distance is that it may be infinite sometimes, if there are too many
%changes between certain pairs of nucleotides. This can be particularly
%noticeable with distances computed from bootstrapped sequences.
%
% Syntax: [D]=dn_logdet(aln)
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
D = zeros(n,n);

for i=1:n-1
for j=i+1:n
	D(i,j) = d_logdet(S(i,:), S(j,:));
	D(j,i) = D(i,j);
end
end


%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function d=d_logdet(seq1, seq2)
 [S,gap] = countntchange(seq1, seq2);
% S=count_ntchanges(seq1, seq2);
 if (det(S)<=0)
	d=inf;
 else
	f1=sum(S);
	f2=sum(S');
	d=(-1/4)*(log(det(S))-(1/2)*log( prod(f1)*prod(f2)));
 end