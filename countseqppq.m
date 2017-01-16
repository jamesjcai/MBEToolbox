function [p1,p2,q]=countseqppq(S)
%COUNTSEQPPQ - Counts transitions (P1 and P2) and transversion (Q)

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (isstruct(S)), S=S.seq; end
[n,m] = size(S);
if (n~=2) error('Not a sequence pair'); end

TS1= [0, 0, 0, 0;
      0, 0, 0, 1;
      0, 0, 0, 0;
      0, 1, 0, 0];

TS2= [0, 0, 1, 0;
      0, 0, 0, 0;
      1, 0, 0, 0;
      0, 0, 0, 0];


TV = [0, 1, 0, 1;
      1, 0, 1, 0;
      0, 1, 0, 1;
      1, 0, 1, 0];

[S] = countntchange(S(1,:), S(2,:));
p1 = sum(sum(TS1.*S));
p2 = sum(sum(TS2.*S));
q = sum(sum(TV.*S));