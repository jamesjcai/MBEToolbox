function [para,lnL]=optimbrchlen2a(treestr,sq,md)

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



options = optimset('fminbnd');
%options = optimset(@fminunc);
%options = optimset('LargeScale','off');
%options=optimset(options,'Display','iter', 'TolX', 1.0000e-008);
LB=0.001;
UB=2.0;
nb = length(brchlen);
f_opt0=inf;

i=0;
while (i<5),
        order=randperm(nb);
	for (k=1:nb),
	        k=order(k);
		[xy,f_opt]=fminbnd(@i_likelifunquick,LB,UB,options,k,brchlen,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);
		%[xy,f_opt]=fminunc(@i_likelifunquick,LB,UB,options,k,brchlen,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);

		
		brchlen(k)=xy(1);
		tolx=f_opt0-f_opt;
		if (tolx>eps), f_opt0=f_opt; end
		if (tolx<eps), i=i+1; end
		disp(sprintf('%5.3f', f_opt0))
	end
end

para=brchlen;
lnL=-f_opt;

disp('brchlen=[0.0356 0.0718 0.0764 0.0861 0.1426 0.0332 0.0280 0.0014]')




function [lnL] = i_likelifunquick(x,k,brchlen,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq)
blen=brchlen;
blen(k)=x;
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



function [lnL] = i_likelifun(x,k,brchlen,treetop,numnode,sq,md)
blen=brchlen;
blen(k)=x;
tr.treetop=treetop;
tr.numnode=numnode;
tr.brchlen=[blen; 0];
[lnL] = -1*treelike(tr,sq,md);