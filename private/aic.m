function [aic] = aic(lnL,k)
%AIC - Computes Akaike's Information Criterion(AIC) from a model
%(AIC) correction (Akaike 1974)
% lnL   - log-likelihood
% k     - number of inferred parameters

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

aic=lnL-k;


function [aicc] = aicc(lnL,k,n)
%AICC - Computes Second-order Akaike (AICC) correction
%(Hurvich and Tsai 1989)
%
% lnL   - log-likelihood
% k     - number of inferred parameters
% n     - sample size

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

aicc=lnL-k-(k(k+1))/(n-k-1);


function [bic] = bic(lnL,k,n)
%BIC - Computes BIC correction (Schwarz 1978)
%
% lnL   - log-likelihood
% k     - number of inferred parameters
% n     - sample size 

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

bic=lnL-(k/2.0)*log(n);