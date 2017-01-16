function [S,N] = estsynnonsynsites(s1,s2);
%GETSYNNONSYNSITES - Estimate Syn- Nonsyn- sites

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[dS,dN,dN_dS,lnL,value] = dc_gy94([s1;s2],1,2);
S=value.S;
N=value.N;
