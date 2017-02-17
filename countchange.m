function [D,gap] = countchange(s1,s2,nword)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if ~(ismember(nword, [4 20 61])),
	error('Wrong NWORD')
end
if (length(s1)~=length(s2)),
	error('Sequences are not of same length.')
end

D=zeros(nword);
WORD=1:nword;

for j=1:nword
   s1sites=(s1==WORD(j));
   for i=1:nword
      D(i,j)=sum(s1sites & (s2==WORD(i)));
   end
end

if (nargout>1),
	gap = length(s1)-sum(D(:));	% gaps
end

% D=accumarray([s1(:),s2(:)],1);    Alex Ellen 2017 homework contribution