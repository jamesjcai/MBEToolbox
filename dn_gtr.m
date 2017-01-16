function [d] = dn_gtr(s1,s2,rmatrix,freq)
%DN_GTR - Distance of GTR model
%
% Syntax: [D,VarD]=dn_gtr(s1,s2,rmatrix,freq)
%
% Inputs:
%    aln    - Alignment structure
%    freq   - (optional) 1x4 vector of equilibrium base frequencies
%
% Outputs:
%    d      - Distance matrix
%
%REF: Yan Z and Nielsen R (2000) Estimating synonymous and nonsynonymous
%     subsititution rates under realistic evolutionary models. p. 34
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[X]=countntchange(s1,s2);

model=modelgtr(rmatrix,freq);

F=((sum(X(:))-trace(X))*model.R)./4;
F=eye(4)*trace(X)./4+F;

PI=diag(model.freq);
%d=-trace(PI*logm(inv(PI)*F));
d=-trace(PI*logm(PI\F));

