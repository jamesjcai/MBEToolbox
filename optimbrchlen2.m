function [para,lnL]=optimbrchlen2(treestr,sq,md)

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


options = optimset('fminbnd');
options=optimset(options,'Display','iter');
LB=0.001;
UB=2.0;


for (k=1:length(brchlen)),
	brchlen
	[xy,f_opt]=fminbnd(@i_likelifun,LB,UB,options,k,brchlen,treetop,numnode,sq,md);
	xy
	brchlen(k)=xy(1);
	brchlen
	pause
end

para=brchlen;
lnL=-f_opt;


function [lnL] = i_likelifun(x,k,brchlen,treetop,numnode,sq,md)
blen=brchlen;
blen(k)=x;
tr.treetop=treetop;
tr.numnode=numnode;
tr.brchlen=[blen; 0];
[lnL] = -1*treelike(tr,sq,md);