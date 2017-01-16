function [p] = stephens85(seq,ptnv)
%STEPHENS85 - Detection of intragenic recombination or gene conversion
%
%
% seq   - sequences
% ptnv  - partition indicator vector
%
%REF: Mol Biol Evol. 1985 Nov;2(6):539-56.
%     Statistical methods of DNA sequence analysis: detection of intragenic 
%     recombination or gene conversion. Stephens JC.

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (nargin<1),    % example input
	seq=randseq(5,100);
end
if (isstruct(seq)), seq=seq.seq; end

[n,m]=size(seq);
if(n<3), error('STEPHENS85:TooFewSeq','Too few sequences.'); end
if(m<3), error('STEPHENS85:TooShortSeq','Sequences are too short.'); end

if (nargin<2),    % example ptnv
	ptnv=[1,1,zeros(1,n-2)];
end

[n2,m2]=size(ptnv);

if (m2~=n), error('STEPHENS85:UnequalLen','Length of partition indicator vector must equal to n.'); end

% there are 2^(m-1)-1 different ways fo dividing m sequences into groups
% nptn=2^(m-1)-1;

p=zeros(1,n2);
for k=1:n2
	inputv = i_seqpartition(seq,ptnv(k,:));
	p(1,k) = i_equation5(inputv);
end


if (nargout<1), 
	disp(' ')
	disp('==========================================');
	disp('Stephens (85) Method: Detection of'); 
	disp('Recombination or Gene Conversion');
	disp('==========================================');
        fprintf ('Number of input sequences (n): %d\n', n);
        fprintf ('Number of ways to dividing them into groups, (2^(n-1)-1): %d\n', 2^(n-1)-1);
        fprintf ('Number of tested phylogenetic partition: %d\n', n2);
	for (k=1:n2),
	fprintf('Phylogenetic partition %d of %d:\n', k, n2);
	ptnv(k,:)
	fprintf('P: %f\n', p(1,k));
	end
	disp('==========================================');
	disp(' ')
% there are 2^(m-1)-1 different ways fo dividing m sequences into groups
% nptn=2^(m-1)-1;
end


function [p] = i_equation5(inputv)
% The following probability is useful for testing whether or not a set of
% s specified sites is clustered relative to the region flanking them.
	p=1;
	n=length(inputv); % length of 
	s=sum(inputv==1);  % number of sites 
	if (s==2)
	      p=(n+1)/3;
	elseif (s>2)
		x=find(inputv==1);
		d0=max(x)-min(x);
		x=1;
		for j=0:(s-2)
		      x=x*((d0-j)/(n-j-1));
		end
		p=(s-(s-1)*(d0+1)/n)*x;               % equation (5)
	else
		p=99.0;
	end
	


function [v] = i_seqpartition(seq,ptnv)
% Returns a vector containing indicators of sites that are compatible 
% to the partition 

%     1     4     1     4     3     4     2     2     1     3
%     2     3     1     4     4     1     2     3     1     2
%     1     2     3     3     3     1     4     3     4     2
%     3     1     1     4     3     2     2     4     3     4
%     2     3     3     1     4     2     4     1     4     4

%ptnv=[1 0 1 1 1];

[n,m]=size(seq);
v=zeros(1,m);

x=find(ptnv==1);
y=find(ptnv==0);

for k=1:m
      site=seq(:,k);
      xx=unique(site(x));
      yy=unique(site(y));
      if (length(xx)==1 && length(yy)==1)
	if (xx~=yy), v(1,k)=1; end
      end      
end

