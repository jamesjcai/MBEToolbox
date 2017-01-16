function [lnL,siteL] = treelike_quick(tr,sq,md,rt,issite)
%TREELIKE - Estimates log-likelihood of a tree
%
% Syntax: [lnL] = likelitree(tree,S,model)
%
% Inputs:
%    tree    - Tree string
%    S       - Aligned sequences
%    model   - Substitution model
%
% Outputs:
%    lnL   - Log-likelihood
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<5), issite=0; end
if (nargin<4), rt=1; end

% the following swap, ensures: treelike(tr,sq,md)==treelike(sq,tr,md)
if (isstr(sq)), temp=tr; tr=sq; sq=temp; end
if (isstruct(sq)), sq=sq.seq; end


[Q] = composeQ(md.R, md.freq);
[W,L] = eig(Q); V = inv(W);
freq = md.freq;
nalpha = length(md.R);

if (isstruct(tr)),
	%treestr=tr.treestr;
	treetop=tr.treetop;
	brchlen=tr.brchlen;
	numnode=tr.numnode;
	%namnode=tr.namnode;
else
	[treetop,numnode,brchlen,namnode]=parsetree(tr);
end


%tr.treestr=treestr;
%tr.treetop=treetop;
%tr.brchlen=brchlen;
%tr.numnode=numnode;
%tr.namnode=namnode;

%[treetop, numnode, brchlen, namnode] = parsetree(tree);
%if ~(n_otu==n_otu2)
%	error('Inconsistent tree and alignment. Need an rooted tree.');
%else
%	inlineL=mbelfcreator(tree,n_otu);
%end

inlineL=mbelfcreator(treetop,numnode);
numbrch = 2*numnode-2;

if ~(length(brchlen)-1==numbrch),
	error('Inconsistent tree and alignment. Need an rooted tree.');
else
	P = zeros(numbrch*nalpha,nalpha);
	for j=1:numbrch
	    P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*brchlen(j)*rt)*V;
	end
end
Lhi = inline(inlineL,'P','numnode','freq','s','nalpha');



if (issite)
	siteL=i_safelog(Lhi(P,numnode,freq,st,nalpha))*sq;
else
	[n_otu2, n_site]=size(sq);
	[patt, npatt, scate] = sitepattern(sq);
	nop=length(npatt);
	lpatt=zeros(1,nop);
	for (k=1:nop),
	    st = patt(:,k);
	    % disp(sprintf('site=%d',k))
	    lpatt(1,k)= i_safelog(Lhi(P,numnode,freq,st,nalpha));
	end
	siteL=lpatt(scate);
end
lnL=sum(siteL);