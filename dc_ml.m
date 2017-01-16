function [dS,dN,dN_dS,lnL] = dc_ml(aln,a,b)
%DC_ML - dS, dN estimation by codeml method
%
% [dS,dN] = dc_ml(aln,a,b)
% calculates dS and dN between sequence a and b in aln.
%
%%

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if nargin < 3,
    a=1;
    b=2;
end

if hasgap(aln)
    aln=rmgaps(aln);
    %aln=rmcodongaps(aln);
end
[dS,dN,dN_dS,lnL] = dc_gy94(aln,a,b);

