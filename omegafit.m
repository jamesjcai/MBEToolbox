function [omg,lnL]=omegafit(sq,tr,kappa)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(sq);
issite=1;
if (m~=1), issite=0; end

%omg=zeros(1,m);
%siteL=zeros(1,m);

[treetop,numnode,brchlen,namnode]=parsetree(tr);
inlineL=mbelfcreator(treetop,numnode);
numbrch = 2*numnode-3;
if ~(length(brchlen)-2==numbrch),
	error('Inconsistent tree and alignment. Need an rooted tree.')
end



options = optimset('fminbnd');
[xy,f_opt]=fminbnd(@i_sitelikelifun,eps,50000000,options,sq,kappa,numnode,inlineL,numbrch,brchlen,issite);
omg=xy(1);
lnL=-f_opt;



function [siteL] = i_sitelikelifun(x,sq,kappa,numnode,inlineL,numbrch,brchlen,issite)
     omega=[x(1)];
     md= modelgy94(omega,kappa);
     [Q] = composeQ(md.R, md.freq);
     [W,L] = eig(Q); V = inv(W);
     freq = md.freq;
     nalpha = length(md.R);
     rate=1;
     siteL=0;

     Lhi = inline(inlineL,'P','numnode','freq','s','nalpha');
     P = zeros(numbrch*nalpha,nalpha);
     for j=1:numbrch
	P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*brchlen(j)*rate)*V; 
     end

     if (issite)
	     st=sq;
	     siteL=-1*log(Lhi(P,numnode,freq,st,nalpha));
     else
	     [n,m]=size(sq);
	     for (k=1:m),
	           st=sq(:,k);
		   siteL=siteL+log(Lhi(P,numnode,freq,st,nalpha));
	     end
  		   siteL=-1*siteL;

     end
