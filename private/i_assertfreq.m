function freq=i_assertfreq(freq)
%I_ASSERTFREQ - Internal function for asserting freq

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

	errortag=0;
	[n,m] = size(freq);
	if ~(n==1&m>=3&m<=4), 
		errortag=1;
	else
	 if (m==4)
		 if (abs(sum(freq)-1)>eps)
			errortag=1;	 
		 end	 	
 	 else
		freq=[freq, 1-sum(freq)];
	 end
	 
	end
	if (errortag==1)
		error('freq supplied conatins error! MBEToolbox uses equal base freq instead.');
		freq=[.25 .25 .25 .25];
	end	
