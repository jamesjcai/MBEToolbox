function [para,lnL]=optimbrchlen1a(treestr,sq,md)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[treetop, numnode, brchlen, namnode] = parsetree(treestr);
%tr.treestr=treestr;
%tr.treetop=treetop;
%tr.numnode=numnode;
%tr.namnode=namnode;

brchlen=ones(numnode*2-2,1)*0.7;	% arbitary starting points
%brchlen=[0.0356 0.0718 0.0764 0.0861 0.1426 0.0332 0.0280 0.0014];





inlineL=mbelfcreator(treetop,numnode);
numbrch = 2*numnode-2;
if (length(brchlen)~=numbrch),
	error('Inconsistent tree and alignment. Need an rooted tree.')
end
Lhi = inline(inlineL,'P','numnode','freq','s','nalpha');
[Q] = composeQ(md.R, md.freq);
[W,L] = eig(Q); V = inv(W);
freq = md.freq;
nalpha = length(md.R);
[n_otu2, n_site]=size(sq);
[patt, npatt, scate] = sitepattern(sq);






options = optimset('fminsearch');
options=optimset(options,'Display','iter');

usebnd=0;

if (usebnd)
	LB=ones(numnode*2-2,1)*0.0005;
	UB=ones(numnode*2-2,1)*1;
	[xy,f_opt]=fminsearchbnd(@i_likelifunquick,brchlen,LB,UB,options,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);
else
	[xy,f_opt]=fminsearch(@i_likelifunquick,brchlen,options,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);
end
para.brchlen=xy;
lnL=-f_opt;



function [lnL] = i_likelifun(x,treetop,numnode,sq,md)
tr.treetop=treetop;
tr.numnode=numnode;
tr.brchlen=[x; 0];
     if (any(x<=eps)|any(x>=1))
    	lnL=inf; return;
     end
     [lnL] = -1*treelike(tr,sq,md);





function [lnL] = i_likelifunquick(x,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq)
%blen=brchlen;
%blen(k)=x;
if (any(x<=eps)|any(x>=1))
    	lnL=inf; return;
end
blen=[x; 0];
rt=1;
	P = zeros(numbrch*nalpha,nalpha);
	for j=1:numbrch
	    P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*blen(j)*rt)*V; 
	end
	nop=length(npatt);
	lpatt=zeros(1,nop);
	for (k=1:nop),
	    st = patt(:,k);
	    xLhi=Lhi(P,numnode,freq,st,nalpha);
		if(xLhi<1e-250)             % too small p
		      lpatt(1,k)=-500;
		else
		      lpatt(1,k)=log(xLhi);
		end
	end
	siteL=lpatt(scate);
lnL=-1*sum(siteL);
