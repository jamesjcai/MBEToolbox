function [lnL] = likelitree(tree,S,model,rate)
%LIKELITREE - Estimates log-likelihood of a tree
%
% Syntax: [lnL] = likelitree(tree,S,model)
%
% Inputs:
% tree - Tree string
% S - Aligned sequences
% model - Substitution model
%
% Outputs:
% lnL - Log-likelihood
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


if (nargin<4), rate=1; end

[Q] = composeQ(model.R, model.freq);
[W,L] = eig(Q); V = inv(W);
freq = model.freq;

%nalpha = 4;
nalpha = length(model.R);

[tree, n_otu, t, nm] = parsetree(tree);
[n_otu2, n_site]=size(S);

if ~(n_otu==n_otu2)
 error('Inconsistent tree and alignment. Need an rooted tree.');
else
 inlineL=mbelfcreator(tree,n_otu);
end

% t = ones(1,n_brc);


n_brc = 2*n_otu - 3;
if ~(length(t)-2==n_brc),
 error('Inconsistent tree and alignment. Need an rooted tree.');
else
 P = zeros(n_brc*nalpha,nalpha);
 for j=1:n_brc
  P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*t(j)*rate)*V;
 end
end

lnL=0;
for (k=1:n_site),
 s = S(:,k);
 Lhi = inline(inlineL,'P','n_otu','freq','s','nalpha');
 lnL = lnL+log(Lhi(P,n_otu,freq,s,nalpha));
end
