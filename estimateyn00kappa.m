function [kappa] = estimateyn00kappa(S)
%ESTIMATEYN00KAPPA - Estimates transition/transversion rate ratio (kappa) by method of YN00
%This function calculates mutational transition/transversion rate ratio
%kappa using 4-fold degenerate sites from pairwise comparisons
%under HKY85, weighting estimates by the numbers of sites
%
%REF: Yan Z and Nielsen R (2000) Estimating synonymous and nonsynonymous
%subsititution rates under realistic evolutionary models. p. 34
%
%Yang and Nielsen (2000) use the fourfold-degenerate sites at the third
%codon positions and the nondegenerate sites to estimate
%$\kappa$. Mutations at the fourfold-degenerate sites do not
%change the amino acid, and thus the transition/transversion
%rate bias at those sites should reflect the mutational
%bias. Mutations at nondegenerate sites all lead to amino
%acid changes and can also be used to estimate
%
% Syntax: [kappa] = estimateyn00kappa(S,method)
%
% Inputs:
%    S        - Sequence pair
%    method   - k2p|f84|hky
%
% Outputs:
%    kappa   - $\kappa=\alpha/\beta$
%
% See also: DN_PDIST

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[n,m]=size(S);
kappa=0;
for (i=1:n),
for (j=i:n),
	if ~(i==j)
	spair =S([i j],:);
	% [t,k] = i_SeqPairDistanceF84(spair)
	[picker0,picker2,picker4] = getDegenerateSitesIndex(spair);
	seq0 = spair(:,picker0);
	seq4 = spair(:,picker4);
	[k0] = i_SeqPairDistanceF84(seq0);
	[k4] = i_SeqPairDistanceF84(seq4);

	wk0=sum(sum(picker0)); wk4=sum(sum(picker4));
	k=(k0*wk0+k4*wk4)/(wk0+wk4);
	kappa=kappa+k/(n*(n-1)/2);
	end
end
end



function [k,t,SEt] = i_SeqPairDistanceF84(S)
%   This calculates kappa and d from P (proportion of transitions) & Q
%   (proportion of transversions) & pi under F84.
%   When F84 fails, we try to use K80.  When K80 fails, we try
%   to use JC69.  When JC69 fails, we set distance t to maxt.
%   Variance formula under F84 is from Tateno et al. (1994), and briefly
%   checked against simulated data sets.

t=99; k=999; SEt=0;

[n,m]=size(S);
[P,Q]=countseqpq(S); P=P./m; Q=Q./m;

% Qsmall=min2(1e-10,0.1/m)

N=[sum(sum(S==1)),sum(sum(S==2)),sum(sum(S==3)),sum(sum(S==4))];
Pi=N./sum(N);		% Pi(A,C,G,T);

tc=Pi(4)*Pi(2);
ag=Pi(1)*Pi(3);

%   Y=p[0]+p[1];    R=p[2]+p[3];  tc=p[0]*p[1]; ag=p[2]*p[3];

Y=Pi(1)+Pi(3);   % Y=TC
R=Pi(2)+Pi(4);   % R=AG

A=tc/Y+ag/R; B=tc+ag; C=Y*R;

a=(2*B+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*C)) - P) / (2*A);
b=1-Q/(2*C);

if (a<=0 || b<=0)
  failF84=1;
  return;
end

a=-.5*log(a); b=-.5*log(b);
if(b<=0)
  failF84=1;
  return;
end

k_F84=a/b-1;
t = 4*b*(tc*(1+ k_F84/Y)+ag*(1+ k_F84/R)+C);
k_HKY = (B + (tc/Y+ag/R)* k_F84)/B;			% k_F84=>k_HKY85

if (k_HKY>999)
	k_HKY=999;
end
k=k_HKY;

if (nargout>2)
     a = A*C/(A*C-C*P/2-(A-B)*Q/2);
     b = A*(A-B)/(A*C-C*P/2-(A-B)*Q/2) - (A-B-C)/(C-Q/2);
     SEt = sqrt((a*a*P+b*b*Q-(a*P+b*Q)^2)/n);
end





function [picker0,picker2,picker4] = getDegenerateSitesIndex(S)

[n,m]=size(S);
S2=codonise64(S);
marker=zeros(n,m);

% deg. tables = {'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA'
% 'ATC' 'ATG' 'ATT' 'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA'
% 'CTC' 'CTG' 'CTT' 'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA'
% 'GTC' 'GTG' 'GTT' 'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA'
% 'TTC' 'TTG' 'TTT'};

degtable=[0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 2 0 2; 0 0 2;
2 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 0; 0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 2;
0 0 4; 0 0 4; 0 0 4; 0 0 4; 2 0 4; 0 0 4; 2 0 4; 0 0 4; 2 0 4; 0 0 4;
2 0 4; 0 0 4; 0 0 2; 0 0 2; 0 0 2; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4;
0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 0; 0 0 2;
0 0 0; 0 0 2; 0 0 4; 0 0 4; 0 0 4; 0 0 4; 0 0 0; 0 0 2; 0 0 0; 0 0 2;
2 0 2; 0 0 2; 2 0 2; 0 0 2; 0 0 0];


for (p=1:n),
for (q=1:m/3)
	deg = degtable(S2(p,q),:);
	a=q;
	q=(q-1)*3+1;	% jump index p every next three bases.
	marker(p,q)=deg(1,1); marker(p,q+1)=deg(1,2); marker(p,q+2)=deg(1,3);
end
end

validflag=zeros(1,m);
k=1:m;
validflag(1,k)=min(marker(:,k))==max(marker(:,k));	% valid columns's flag, 0 or 1

pick0=marker==0;
pick2=marker==2;
pick4=marker==4;
picker0=pick0(1,:)&validflag;
picker2=pick2(1,:)&validflag;
picker4=pick4(1,:)&validflag;



%function [P,Q]=countseqpq(S)
%
%[n,m] = size(S);
%if (n~=2) error('Not a sequence pair'); end
%
%TS = [0, 0, 1, 0;
%      0, 0, 0, 1;
%      1, 0, 0, 0;
%      0, 1, 0, 0];
%
%TV = [0, 1, 0, 1;
%      1, 0, 1, 0;
%      0, 1, 0, 1;
%      1, 0, 1, 0];
%
%P = zeros(n);
%Q = zeros(n);
%	[S,gap] = countntchange(S(1,:), S(2,:));
%	P = sum(sum(TS.*S));
%	Q = sum(sum(TV.*S));