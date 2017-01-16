function [stpchain] = gibbsmotif(s,w,nm)
%GIBBSMOTIF - A Very Basic Gibbs Sampler for Motif Detection
%
% s - sequences
% w - motif width
% nm - how many motif wanted

%Lawrence, C.E., Altschul, S.F., Boguski, M.S., Liu, J.S., Neuwald, A.F. and 
%Wootton, J.C. 1993. Detecting subtle sequence signals: a Gibbs sampling 
%strategy for multiple alignment. Science 262: 208-214.

%%example file:
%>SeqName1_or_space
%gtacacacacacgtacgtgtagtcagtactacacgtacgt
%gtcagtcgacacgagtacgtactacttcgcgcgtacgtac
%>SeqName2_or_space
%caatacgtaactgacacttgcagtagtcagtactgtgtacgt
%>SeqName3_or_space
%gtacacgtcaactactgtagtagtcagtactacatgactg
%gtatatactcggcgcgtcatactactcgtcgtactacgta
%actgaccagagtacacgacgtacgtacgtacacgtactac

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (nargin<2)
	%s=[randseq(30);randseq(30);randseq(30);randseq(30);randseq(30);...
	%   randseq(30);randseq(30);randseq(30);randseq(30);randseq(30)];
	w=20;
	nm=3;
end
[n,m]=size(s);


%2. Calculate the background frequencies of A,C,G,T from all the sequences
[dump,freq]=ntcomposition(reshape(s,1,n*m));

leftoutid=1;
[score0,stp]=calcscore(s,m,n,w,freq,leftoutid);


nstep=1000;
stpchain=zeros(nstep,n);

while (nstep>0),
	leftoutid=ceil(rand*(n-1));
	[score,stp]=calcscore(s,m,n,w,freq,leftoutid);
	if (score>=score0)		
	      stpchain(nstep,:)=stp;			% jump to new start;
	      nstep=nstep-1;
	else
		if (rand>score/score0),
			stpchain(nstep,:)=stp;		% jump to new start;
		      nstep=nstep-1;
		end
	end
end




%showseqbolck(s,stp,w);


%6. Randomly generate another start position of the motif for that left-out sequence.
%7. Score that sequence with its new start position.
%8. Compare this new score with its original score.
%9. If newscore >= oldscore, then jump to that new start position, else jump to that 
%   new start position with probability = 



function [score,stp] = calcscore(s,m,n,w,freq,id)

%3. Generate random start positions for the motif in each sequence.
stp=ceil(rand(1,n)*(m-w));   % startpoints
x=[stp;stp+w-1];

sb=[];
for (i=1:n),
      sb=[sb;s(i,[x(1,i):x(2,i)])];
end

%4. Construct position specific score matrix from all sequences except one.
leftout=sb(id,:);
sb(id,:)=[];

pssm=zeros(4,w);
for (k=1:4), pssm(k,:)=sum(sb==k,1); end
pssm=pssm./n;

%5. Score the left-out sequence according to the position specific score matrix:
Pi=freq(leftout);
Pki=zeros(1,w);
for (k=1:w),
      col=pssm(:,k);
      Pki(k)=col(leftout(k));
end
score=sum(Pki./Pi);



function showseqbolck(s,stp,w)

x=[stp;stp+w-1];
NT='ACGT';
[n,m]=size(s);
sb=[];
for (i=1:n),
      sb=[sb;s(i,[x(1,i):x(2,i)])];
end
NT(sb)
