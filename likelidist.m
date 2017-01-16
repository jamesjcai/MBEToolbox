function [lnL] = likelidist(t,model,s1,s2)
%LIKELIDIST - Estimates log-likelihood of branch length (distance)
%Estimates log-likelihood of distance/branch length given subst. model
%This function is the same as SEQPAIRLIKELI, just a different implementation
%
% Syntax: [lnL] = likelidist(t,model,s1,s2)
%
% Inputs:
%    t       - Distance/branch length
%    model   - Substitution model
%    s1      - Nucleotide sequence 1
%    s2      - Nucleotide sequence 2
%
% Outputs:
%    lnL   - Log-likelihood
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



R=model.R;
freq=model.freq;

%L=1;
PI=diag(freq);
Q=composeQ(R,PI);

[V,D] = eig(Q*t);
P=V*diag(exp(diag(D)))/V;
%P=expm(Q*t);

[patt,npatt] = sitepattern([s1;s2]);

lnL=0;
for k=1:length(npatt)
      b1=patt(1,k); b2=patt(2,k);
      %L = L * power(P(b1,b2)*freq(b1),npatt(k));
      p=P(b1,b2);
      % if (p<1.0e-50), x=eps; end % too samll p
      lnL = lnL + log(p*freq(b1))*npatt(k);
end
