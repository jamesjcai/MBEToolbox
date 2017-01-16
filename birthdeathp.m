function [p]=birthdeathp(lamda,t,s,c)
%BIRTHDEATHP - transition probability of mirth and death model
%
%
%
%REFS: Hahn MW (2005) and Bailey (1964)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


lamda=0.002;
t=1.2;
s=1;
c=11;

alpha=(lamda.*t)./(1+lamda.*t);
p=0;
for j=0:min([s,c])
    p=p+nchoosek(s,j)*nchoosek(s+c-j-1,s-1)*power(alpha,(s+c-2*j))*power((1-2*alpha),j);
end

