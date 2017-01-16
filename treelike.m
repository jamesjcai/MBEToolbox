function [lnL,siteL] = treelike(tree,seq,model,rate,issite)
%TREELIKE - Estimates log-likelihood of a tree
%
% Syntax: [lnL] = treelike(tree,seq,model)
%
% Inputs:
% tree - Tree string
% seq - Aligned sequences
% model - Substitution model
% rate - Sub. rate
%
% Outputs:
% lnL - Log-likelihood
% siteL - Log-likelihoods of sites
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
if (nargin<4), rate=1; end

% the following swap, ensures: treelike(tree,seq,model)==treelike(seq,tree,model)
if (ischar(seq)), temp=tree; tree=seq; seq=temp; end
if (isstruct(seq)), seq=seq.seq; end


[Q] = composeQ(model.R, model.freq);
[W,L] = eig(Q); V = inv(W);
freq = model.freq;
nalpha = length(model.R);

if (isstruct(tree)),
 %treestr=tree.treestr;
 treetop=tree.treetop;
 brchlen=tree.brchlen;
 numnode=tree.numnode;
 %namnode=tree.namnode;
elseif (ischar(tree)),
 [treetop,numnode,brchlen,namnode]=parsetree(tree);
else
 error('something wrong with tree string')
end


%tree.treestr=treestr;
%tree.treetop=treetop;
%tree.brchlen=brchlen;
%tree.numnode=numnode;
%tree.namnode=namnode;

%[treetop, numnode, brchlen, namnode] = parsetree(tree);
%if ~(n_otu==n_otu2)
% error('Inconsistent tree and alignment. Need an rooted tree.');
%else
% inlineL=mbelfcreator(tree,n_otu);
%end

inlineL=mbelfcreator(treetop,numnode);
numbrch = 2*numnode-2;

if (length(brchlen)-1~=numbrch),
 error('Inconsistent tree and alignment. Need an rooted tree.')
else
 P = zeros(numbrch*nalpha,nalpha);
 for j=1:numbrch
  P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*brchlen(j)*rate)*V;
   %P((j-1)*nalpha+1:j*nalpha,:)=expm(Q*brchlen(j)*rate);
 end
end


%Lhi = inline(inlineL,'P','numnode','freq','s','nalpha');
Lhi = inline(inlineL,'P','freq','s','nalpha');

if (issite)
 st=seq; % should check it valid
 siteL=log(Lhi(P,freq,st,nalpha));
else
 [n_otu2, n_site]=size(seq);
 [patt, npatt, scate] = sitepattern(seq);
 nop=length(npatt);
 lpatt=zeros(1,nop);
 for (k=1:nop),
  st = patt(:,k);
 % disp(sprintf('site=%d',k))
 % Lhi(P,numnode,freq,st,nalpha)
 xLhi=Lhi(P,freq,st,nalpha);
 if(xLhi<1e-250) % too small p
  lpatt(1,k)=-500;
 else
  lpatt(1,k)=log(xLhi);
 end

 if (nargout<1),
  disp(sprintf('site=%4d, no=%4d x lnL=%3.4f',k,npatt(k),lpatt(k)))
 end

 end

 siteL=lpatt(scate);
% lpatt;
% sum(siteL);
end
lnL=sum(siteL);
