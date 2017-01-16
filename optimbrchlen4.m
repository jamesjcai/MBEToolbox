function [para,lnL]=optimbrchlen4(treestr,sq,md)
% The Gibbs sampler only changes one random variable at a time
% - Slow convergence
% - High-probability states may not be reached because reaching them requires 
%   going through low-probability states

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


[lnL0] = i_likelifunquick(brchlen(1),1,brchlen,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);


%[lnL0] = i_likelifun(brchlen(1),1,brchlen,treetop,numnode,sq,md);
omega=.1;
nstep=5000;

for (j=1:nstep),
%disp(int2str(j))
for (k=1:length(brchlen)),		
		newlen=brchlen(k)-omega/2+omega*rand;
		if (newlen<0), newlen=-1*newlen; end
		[lnL] = i_likelifunquick(newlen,k,brchlen,numbrch,nalpha,W,L,V,patt,npatt,scate,Lhi,numnode,freq);
		alpha=min(1, lnL/lnL0);
		if (rand<alpha), brchlen(k)=newlen; end
end
end




para=brchlen;
lnL=i_likelifun(brchlen(1),1,brchlen,treetop,numnode,sq,md);




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