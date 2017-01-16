function [pssm,freq] = buildpssm(s)
%BUILDPSSM - Builds a position specific scoring matrix (PSSM)
%Position specific scoring matrix (PSSM) or Positional weight matrix (PWM)
%is a motif descriptor. It attempts to capture the intrinsic variability
%characteristic of sequence patterns. A Profile it is usually derived from a
%set of aligned sequences functionally related. For instance, below we have the
%sequence of ten vertebrate donor sites, aligned at the boundary exon/intron.
%
% Syntax: [pssm,freq,LnF] = buildpwm(s)
%
% Inputs:
%    s     - Sequences aligned
%
% Outputs:
%    pssm     - Raw numbers for A C G T
%    freq     - Frequencies
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



[n,m]=size(s);

pssm=zeros(4,m);
for (k=1:4),
      pssm(k,:)=sum(s==k,1);
end

freq=pssm./n;