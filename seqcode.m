function [NT,AA] = seqcode()
%SEQCODE - Return vector for mapping sequence letters to integers
%
% Syntax: [NT,AA] = seqcode
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


NT = 'ACGT-';
if (nargout>1)
AA = 'ARNDCQEGHILKMFPSTWYV*-';
end

% AANames = {'ala' 'arg' 'asn' 'asp' 'cys' 'gln' 'glu' 'gly' 'his' 'ile' 'leu' 'lys' 'met' ...
%            'phe' 'pro' 'ser' 'thr' 'trp' 'tyr' 'val'};