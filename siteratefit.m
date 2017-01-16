function [r,lnL]=siteratefit(s,tr,md)
%SITERATEFIT - Site-specific rate estimate for a given tree and site. 

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(s);
if (m~=1), error('Please provide single site.'); end
[r,lnL]=optimsiterate(s,tr,md);


function [r,lnL]=optimsiterate(s,tr,md)
options = optimset('fminsearch');
[xy,f_opt]=fminsearch(@negloglike,5.0,options,s,tr,md);
r=xy(1);
lnL=-f_opt;


function [nll] = negloglike(x,s,tr,md)
    rate=[x(1)];
    if (x(1)>=20 | x(1)<eps)
    	nll=Inf; return;
    else
         nll = -1*siteratelike(s,tr,md,rate);   
    end
     
