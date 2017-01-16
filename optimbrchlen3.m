function [para,lnL]=optimbrchlen3(treestr,sq,md)
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


[lnL0] = i_likelifun(brchlen(1),1,brchlen,treetop,numnode,sq,md);
omega=.1;
nstep=5000;

for (j=1:nstep),
%disp(int2str(j))
for (k=1:length(brchlen)),		
		newlen=brchlen(k)-omega/2+omega*rand;
		[lnL] = i_likelifun(newlen,k,brchlen,treetop,numnode,sq,md);
		alpha=min(1, lnL/lnL0);
		if (rand<alpha), brchlen(k)=newlen; end
end
end

para=brchlen;
lnL=i_likelifun(brchlen(1),1,brchlen,treetop,numnode,sq,md);


function [lnL] = i_likelifun(x,k,brchlen,treetop,numnode,sq,md)
blen=brchlen;
blen(k)=x;
tr.treetop=treetop;
tr.numnode=numnode;
tr.brchlen=[blen; 0];
[lnL] = -1*treelike(tr,sq,md);