function [model] = modelhky(freq,kappa)
%MODELHKY - Returns a structure of model HKY
%\paragraph{The HKY model}
%The Hasegawa, Kishino and Yano (HKY) model (Hasegawa et al. 1985) allows for
%a different rate of transitions and transversions as well as unequal frequencies
%of the four nucleotides (base frequencies). The parameters requires by this model
%are transition to transversion ratio (TS/TV) and the base frequencies. There are
%a number of simpler models that are specific cases of the HKY modesl. If the base
%frequencies are set equal then the model becomes equivalent to the Kimura 2-parameter
%(K2P) model (Kimura, 1980). If the TS/TV is set to 0.5 as well, then it becomes
%equivalent to the Jukes-Cantor (JC69) model (Jukes and Cantor, 1969). If the TS/TV is
%et to 0.5 and the base frequencies are not equal then the model is equivalent to the
%F81 model (Felsenstein, 1981).
%
% Syntax: [model] = modelhky(freq,kappa)
%
% Inputs:
%    kappa   - $\kappa=\alpha /\beta$
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


if (nargin<2)
	warning('kappa does not supplied. Assuming kappa=1')
    kappa=1;
end

if (nargin<1)
	warning('freq does not supplied. Assuming equal freq')
    freq=[.25 .25 .25 .25];
end


%    Transition/transversion ratio = 2 (kappa = 4)
%    Assumed nucleotide frequencies (set by user):

i_assertfreq(freq);
R=i_hky_ttr(kappa);

a=R./(R+1);
b=0.5*(1./(R+1));
rmatrix=[b a b b a b];
model=model_nt('hky','rmatrix',rmatrix,'freq',freq);

function [R] = i_hky_ttr(kappa)
R=kappa./2;