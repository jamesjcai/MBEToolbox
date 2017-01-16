function [pos] = randcut(G,N,len)
%RANDCUT - Cuts genome DNA randomly into shot fragment
%Genome of size, G, is cutted random into N fragments of average length, len.

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (len>=G), error('No cut.'); end

pos=zeros(N,2);
for (k=1:N),
      ps=1+floor(rand*G);
      pe=ps+len+floor(randn*100);
      pe=min(pe,G);
      pos(k,1)=ps;
      pos(k,2)=pe;
      %fprintf('%d - %d\n',ps,pe);
end