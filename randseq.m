function [seq,seqstr]=randseq(n,m,seqtype)
%RANDSEQ - Generates a random sequence
%
% Syntax: [seq]=randseq(n,m)
%
% Inputs:
%    m         - Length of seq
%    n         - Number of seq
%    seqtype   - DNA=1,Coding DNA=2,Protein=3
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


if (nargin<3), seqtype=1; end
if (nargin<2), m=30; end
if (nargin<1), n=3; end

switch (seqtype)
    case (1)
         nalpha=4;
	 calpha='ACGT';
    case (2)
         nalpha=61;
 	 calpha='ACGT';
    case (3)
         nalpha=20;
	 calpha='ARNDCQEGHILKMFPSTWYV';
end



if (seqtype==2), m=round(m/3); end;
seq=ceil(rand(n,m)*nalpha);
if (seqtype==2), [seq] = decodonise61(seq); end;



if (nargout>1),
	seqstr=calpha(seq);
end



function [seq,nt]=randseq_old(n,freq)
	if nargin<2
	    freq=[.25, .25, .25, .25];
	end
	if nargin<1
	    n=100;
	    disp ('100 base generated.')
	end

	if sum(freq)~=1
	    freq(4)=1-sum(freq([1:3]));
	end

	counter = fix(n*freq);
	seq=[];

	for k=1:4
	    seq=cat(2,seq,ones(1,counter(k))*k);
	end

	for k=1:n-length(seq)
	   seq=cat(2,seq,fix(rand*4+1));
	end
	seq=seq(randperm(length(seq)));

	if (nargout>1),
	      NT='ACGT';
	      nt=NT(seq);
	end
