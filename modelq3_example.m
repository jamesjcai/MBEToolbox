function [para,lnL]=modelq3_example

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

seq=[3     2     3     2
     2     1     2     2
     2     1     1     2];        % example data
 
tree='((a:0.1,b:0.2):0.4,c:0.1);'; % fixed tree;
freq=[0.3,0.4,0.3];    % equibrium frequencies


disp('Before optimization:')
r=[1 2 3 3 2 1]
model=modelq3(r,freq);
[lnL,siteL] = treelike(tree,seq,model)


[para,lnL]=i_optimt(seq,tree,freq);
disp('After optimization:')
r=[para,1]
lnL




function [para,lnL]=i_optimt(seq,tree,freq)
	er=[1 1 1 1 1];  % initial values for r
	options = optimset('fminsearch');	
	[para,f_opt]=fminsearch(@i_likelifunt,[er],options,seq,tree,freq);
	lnL=f_opt;

    
function [lnL] = i_likelifunt(x,seq,tree,freq)
	lnL=inf;
    if any(x)<0, return; end
	a=x(1); b=x(2); c=x(3); d=x(4); e=x(5);
    r=[a b c d e 1];
    model=modelq3(r,freq);
    [lnL] = treelike(tree,seq,model);
    
    

    
    