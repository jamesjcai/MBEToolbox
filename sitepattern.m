function [patt, npatt, scate] = sitepattern(S)
%SITEPATTERN - Infers site patterns for a given alignment
%
% [patt, npatt, scate] = sitepattern(S)
%
%Same as sitecombcrunch in PAML

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[n,m] = size(S);

patt=unique(S','rows');
patt=patt';
p=size(patt,2);
npatt=zeros(1,p);
scate=zeros(1,m);

for (i=1:p),
for (j=1:m),
      if (isequal(patt(:,i), S(:,j)))
          npatt(1,i)=npatt(1,i)+1;
          scate(1,j)=i;
      end
end
end

% assert:  patt(:,scate)==S;