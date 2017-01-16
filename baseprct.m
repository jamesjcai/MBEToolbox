function p=baseprct(seq,dim)
%BASEPRCT - percentage of base in sequences
%
% dim=1, row-wise; dim=2, column-wise

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if nargin<2, dim=1; end

[n,m]=size(seq);
if isstr(seq)
    %sum(seq=='A'|seq=='C'|seq=='G'|seq=='T',2)
    seq=upper(seq);
    if dim==1
        p=sum(ismember(seq,'ATCG'),2)./m;
    else
        p=sum(ismember(seq,'ATCG'),1)./n;
    end
end