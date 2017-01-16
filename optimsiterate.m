function [r,lnL]=optimsiterate(st,tr,md)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(st);
if (m~=1), error('Please provide single site.'); end

%r=zeros(1,m);
%siteL=zeros(1,m);

[treetop,numnode,brchlen,namnode]=parsetree(tr);
inlineL=mbelfcreator(treetop,numnode);
numbrch = 2*numnode-2;

%if ~(length(brchlen)-2==numbrch),
%	error('Inconsistent tree and alignment. Need an rooted tree.')
%end

[Q] = composeQ(md.R, md.freq);
[W,L] = eig(Q); V = inv(W);
freq = md.freq;
nalpha = length(md.R);



options = optimset('fminbnd');
[xy,f_opt]=fminbnd(@i_sitelikelifun_quicker,eps,20,options,st,numnode,freq,inlineL,numbrch,nalpha,brchlen,L,W,V);
%[xy,f_opt]=fminbnd(@i_sitelikelifun_slower,eps,20,options,st,numnode,freq,inlineL,numbrch,nalpha,brchlen,L,W,V,tr,md);
r=xy(1);
lnL=-f_opt;

function [siteL] = i_sitelikelifun_slower(x,st,numnode,freq,inlineL,numbrch,nalpha,brchlen,L,W,V,tree,model)
     rate=[x(1)];
     siteL=log(treelike(tree,st,model,rate));


function [siteL] = i_sitelikelifun_quicker(x,st,numnode,freq,inlineL,numbrch,nalpha,brchlen,L,W,V)
     rate=[x(1)];
     Lhi = inline(inlineL,'P','numnode','freq','s','nalpha');
     P = zeros(numbrch*nalpha,nalpha);
     for j=1:numbrch
	P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*brchlen(j)*rate)*V; 
     end
     x=Lhi(P,numnode,freq,st,nalpha);
     x=max(x,eps);
     siteL=-1*log(x);

