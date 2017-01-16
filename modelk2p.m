function [model] = modelk2p(kappa)
%MODELK2P - Returns a structure of model K2P
%\paragraph{The K2P model}
%The K2P model was described by Kimura (1980). It is like the JC model in
%assuming equal base frequencies, but allows the rate of transition-type
%substitutions to differ from the rate of transversion-type substitutions.
%The ratio of these two instantaneous rates is $\kappa$, and you will be free
%to vary both $\kappa$ and $\mu t$ when this model is selected. Setting
%$\kappa$ equal to 1.0 manually makes this model identical with the JC model.
%The base frequency parameters are inactivated under this model, as they are
%forced to be equal.
%
%Transition-transversion ratio:
%$$R=\kappa/2$$
%
%Evolutionary distance:
%$$d=(\frac{\kappa+2}{4})\mu t$$
%
%\paragraph{$\kappa$ parameter}
%The $\kappa$ parameter represents the ratio of the instantaneous rate of
%transition-type substitutions to transversion-type substitutions. It assumes
%the value 1.0 for models in which all substitutions are taken to occur at the
%same rate (i.e., the JC and F81 models). In the K2P and HKY models, the rate of
%transversions is $\beta$, with the rate of transitions being determined as
%$\alpha=\kappa \beta$.
%
% Syntax: [model] = modelk2p(kappa)
%
% Inputs:
%    kappa   - $\alpha=\kappa \beta$
%
% Outputs:
%    model.R      - Rate matrix
%    morel.freq   - Equilibrium frequency parameters, 1x4
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


if (nargin<1)
	rmatrix=[];
else
	if (isreal(kappa)&kappa>=0&kappa<realmax)
		R=kappa/2;
		a=R./(R+1);
		b=0.5*(1./(R+1));
		rmatrix=[b a b b a b];
	end
end

model=model_nt('k2p','rmatrix',rmatrix);

%    Transition/transversion ratio = 2 (kappa = 4)
%    Assumed nucleotide frequencies (set by user):