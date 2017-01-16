function [para,lnL]=optimbrchlen(treestr,sq,md)

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

options = optimset('fminsearch');
options=optimset(options,'Display','iter');


usebnd=1;

if (usebnd)
	LB=ones(numnode*2-2,1)*0.001;
	UB=ones(numnode*2-2,1)*1;
	[xy,f_opt]=fminsearchbnd(@i_likelifun,brchlen,LB,UB,options,treetop,numnode,sq,md);
else
	[xy,f_opt]=fminsearch(@i_likelifun,brchlen,options,treetop,numnode,sq,md);
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