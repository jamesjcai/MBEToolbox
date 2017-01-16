function [s2]=mutateseq(s1,model,t)
%MUTATESEQ - Monte Carlo simulation of DNA and amino acid sequence evolution
%along phylogenetic trees

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


%rmatrix=[1.0 1.33333 1.0 1.0 1.333333];
%freq=[.1 .2 .3 .4];
%model = modelgtr(rmatrix,freq);
%aln.seq=[s1;s2];
%aln.seqnames={'s1','s2'};
%aln.seqtype=1;
%aln.geneticcode=0;
%dn_jc(aln)
%model=modeljtt;
%s2=mutateseq(s1,model,0.8);

[n,m]=size(s1);
if (n~=1), error('Please supply single sequence'); end

R=model.R;
PI=model.freq;
Q = composeQ(R,PI);

[V,D] = eig(Q*t);
P=V*diag(exp(diag(D)))/V;
%P=expm(Q*t);

state_gap=length(R)+1;
s2=ones(n,m)*state_gap; % fill s2 with gaps(unknown states).
for (k=1:m),
      P_k=P(:,s1(1,k));
      s2(1,k)=SampleFromMultinomial(P_k);
end




function zzz=SampleFromMultinomial(p)
% SampleFromMultinomial - Samples from a multinomial distribution
% If the elements of p do not sum to one, p gets normalised
%
% INPUT
% p(1:K) - array of probabilities
%
% OUTPUT
% zzz - Element in (1,2,...,K) = sample from the multinomial
%
% USAGE
% zzz=SampleFromMultinomial(p)

% -- Error checking --
if sum(p<0) > 0
	error('p contains negative elements');
end
Z=sum(p);
if Z==0
	error('All elements in p are zero.')
end
if Z~=1
	p=p/Z;
end

%K= length(p);
%cumulative_dist=zeros(1,K);
%for i=1:K-1
%	cumulative_dist(i+1)=cumulative_dist(i)+p(i);
%end
%% This is an array with treshold(1)=0, cumulative_dist(K)=1-p(K)

p=p';
cumulative_dist=cumsum([0,p(1:end-1)]);

x=rand;
zzz=sum(x>=cumulative_dist);
