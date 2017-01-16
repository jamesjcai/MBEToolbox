function [seq2] = seqgen(n,m,seqtype)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if nargin<3, seqtype=1; end
if nargin<2, m=300; end
if nargin<1, n=10; end

    
seq1=randseq(1,m,seqtype);
switch seqtype
    case 1
        model=modeljc;
    case 2
        model=modelgy94;
    case 3
        model=modeljtt;
end

seq2=[];
for k=1:n
    t=abs(randn);
    sq=mutateseq(seq1,model,t);
    seq2=[seq2;sq];
end
