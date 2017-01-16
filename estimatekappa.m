function [kappa] = estimatekappa(S,method)
%ESTIMATEKAPPA - Estimates kappa for given sequence pair
%
% Syntax: [kappa] = estimatekappa(S,method)
%
% Inputs:
%    S        - Sequence pair
%    method   - k2p|f84|hky
%
% Outputs:
%    kappa   - $\kappa=\alpha/\beta$
%
% See also: DN_PDIST

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2)
      method='f84';
end

[n,m]=size(S);
if (n~=2) error('Not a sequence pair'); end
if (m<2) error('Too short sequence'); end


largek=999;
[P1,P2,Q]=countseqppq(S);
P1=P1./m; P2=P2./m; Q=Q./m;

switch (lower(method))
    case {'f84','hky'}
	N=[sum(sum(S==1)),sum(sum(S==2)),sum(sum(S==3)),sum(sum(S==4))];
	Pi=N./sum(N);		% Pi(A,C,G,T);
	tc=Pi(4)*Pi(2);
	ag=Pi(1)*Pi(3);
	R=Pi(1)+Pi(3);
	Y=Pi(2)+Pi(4);
end

switch (lower(method))
    case {'k2p'}
      a1=1-2*(P1+P2)-Q;   b=1-2*Q;
      a1=-log(a1);  b=-log(b);
      a1=0.5*a1-0.25*b; b=0.25*b;
      if(b>eps),
	kappa = a1./b;
      else
	kappa=largek;
      end
    case {'f84'}
      a1=(2*(tc+ag)+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*Y*R)) -P1-P2) / (2*tc/Y+2*ag/R);
      b = 1 - Q/(2*Y*R);
      a1=-log(a1); b=-log(b);
      a1=.5*a1;  b=.5*b;

      if(b>eps),
        kappa = a1./b-1;
      else
	kappa=largek;
      end
    case {'hky'}
      kappa = largek;
      a1=1-Y*P1/(2*tc)-Q/(2*Y);
      a2=1-R*P2/(2*ag)-Q/(2*R);
      b=1-Q/(2*Y*R);
      if (a1<=0 || a2<=0 || b<=0),
           return;
      end
      a1=-log(a1); a2=-log(a2); b=-log(b);
      a1 = -R/Y*b + a1/Y;
      a2 = -Y/R*b + a2/R;
      if (b>0), kappa = min((a1+a2)./(2*b), largek); end
end