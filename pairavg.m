function [d] = pairavg(funfcn,aln,varargin)
%PAIRAVG - Pairwise average.  
%   USAGE: [d] = pairavg(funfcn,aln)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

funfcn = fcnchk(funfcn,length(varargin));
%if ~isa(funfcn,'function_handle'), ...

[n,m]=size(seq);
res=nchoosek(n,2);
k=1;
for i=1:n-1
for j=i+1:n
    res(k)=feval(funfcn,seq([i,j],:),varargin{:});
    k=k+1;
end
end

d=mean(res);