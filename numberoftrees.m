function [ntree] = numberoftrees(ns,rooted)
%NUMBEROFTREES - Number of rooted or unrooted trees

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<1), error('Please specific ns.'); end
if (nargin<2), rooted=1; end
if (ns<4), error('ns too small!'); end
if (ns>15), error('ns too large!'); end
ntree=1;

for (i=4:ns),
  ntree=ntree*(2*i-5);
end
if (rooted), ntree=ntree*(2*i-3); end