function [D]=dn_jin_nei90(aln,alpha)
%DN_JIN_NEI90 - Jin-Nei Gamma distance
%This method applies to nucleotides only. This again uses transition and
%transversion rates. The shape parameter, i.e. "alpha", is the square of the
%inverse of the coefficient of variation of the average substitution.
%
% Syntax: [D]=dn_jin_nei90(aln,alpha)
%
% Inputs:
%    aln     - Alignment structure
%    alpha   - Shape parameter of gamma distribution
%
% Outputs:
%    D      - Distance matrix
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


if(nargin==1), alpha=1.0; disp('Using defalt alpha=1.0'); end

[n,m] = size(aln.seq);
[P,Q]=countseqpq(aln);
P=P./m; Q=Q./m;

% ??? l = (p/m)+2*(q/m);

% Jin-Nei Gamma distance
% This method applies to nucleotides only. This again uses transition
% and transversion rates. As with the Kimura two parameter method, gaps
% and ambiguous symbols other than R and Y are not oncluded in the score.
% The shape parameter, i.e. "a", is the square of the inverse of the
% coefficient of variation of the average substitution,
%

% L = average substituition = transition_rate + 2 * transversion_rate
% a = (average L)^2/(variance of L)

% P = transitions/npos
% Q = transversions/npos
% npos - number of positions scored

%d = 0.5 * a * ((1-2*p-q)^(-1/a) + 0.5 * (1-2*q)^(-1/a) -3/2) ;

W1=1-2*P-Q;
W2=1-2*Q;
D = 0.5*alpha*(W1.^(-1/alpha) + 0.5*W2.^(-1/alpha)-3/2) ;

% It is suggested [Jin et al.], in general, that the distance be calculated with an a-value of 1.
% However, the user can specify their own value, using the "-parametera" option, or calculate for
% each pair of sequence, using "-calculatea".
%
% Reference:
% L. Jin and M. Nei, Mol. Biol. Evol. 1990, 7, 82.
%
