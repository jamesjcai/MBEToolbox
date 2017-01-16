function [s] = revseq(s0)
%REVSEQ - Return reverse strand of nucleotide sequences
%
% Syntax: [s] = revseq(s0)
%
% Inputs:
%    s0     - Sequences
%
% Outputs:
%    s    - New sequences
%
% See also: TRANSLATESEQ, REVCOMSEQ

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (isstruct(s0))
      if ~(isvalidaln(s0,'NUCLEOTIDE'))
         error ('ERROR: Not nucleotide sequence.')
      end
         s0=s0.seq;
end
%s=s0(:,end:-1:1);
s = fliplr(s0);
