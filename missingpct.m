function [p]=missingprct(aln,isnuc)
%MISSINGPCT - percentage of missing bases (amino-acids) in sequence
%
% [p]=missingpct(aln,isnuc)
% isnuc = 1 (default)  % sequence type is nucleotide
% isnuc = 0            % sequence type is amino-acid

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if nargin<2
    isnuc=1;
end

if isstruct(aln)
    seq=aln.seq;
else
    seq=aln;
end

if isnuc
idx=aln.seq>0 & aln.seq<5;
p=sum(idx,2)./size(idx,2);
else
idx=aln.seq>0 & aln.seq<20;
p=sum(idx,2)./size(idx,2);
end