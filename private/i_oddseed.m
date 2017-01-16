function [seed]=i_oddseed()
% internal function [seed]=i_oddseed()
% generate odd random number
%
% Matlab MBEToolbox internal function

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

	seed = round(1000*rand);
if (~(mod(seed,2)==1))
      seed=seed+1;
end
