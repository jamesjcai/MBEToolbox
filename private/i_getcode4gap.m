function [n] = i_getcode4gap(seqtype)

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005


switch (upper(seqtype))
    case {'PROTEIN'}
         [~,AA] = seqcode;
	 n=find(AA=='-');
         % n=22;
    case {'DNA'}
         [NT] = seqcode;
	 n=find(NT=='-');
         % n=5;
    otherwise
	error('Wrong seqtype.')
end

