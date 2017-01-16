function [inv]=invdistrib(f,se)
%INVDISTRIB - Invariable sites model
%construct discrete rate distribution with two rates
%(one invariable and one variable) (two-rate model with mean rate = 1.0)
%
% Syntax: [inv]=invdistrib(f)
%
% Inputs:
%    f    - Fraction of invariable sites
%    se   - S.E.
%
% Outputs:
%    inv   - Invariable sites model
%
% See also: GAMMADISTRIB

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


showSE=0;
if (nargin==2)
  showSE=1;
end


% number of rate categories
inv.numRates=2;

% rates of each rate category
inv.rate=[0, 1./(1-f)];         % ensures that mean rate = 1.0, i.e., sum(inv.rate.*inv.prob)=1

% probability of each rate
inv.prob=[f, 1-f];

if (showSE)
	inv.fracSE = se;
end


fprintf('Model of rate heterogeneity: Invariable sites model\n');
fprintf('Number of rate categories: %d\n', inv.numRates);
fprintf('Fraction of invariable sites: %2.2f', f);

if (showSE)
	fprintf(' (S.E. %2.2f)\n', inv.fracSE);
end