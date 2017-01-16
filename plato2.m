function [Q] = plato2(sq,tr,md)
%PLATO - Detecting recombination with window methods
%Partial Likelihoods Assessed Through Optimisation

%Detecting recombination with window methods
%REF: Grassly, Holmes (Mol. Biol. Evol. 14)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (isstruct(sq)), sq=sq.seq; end
%if (nargin<5), sp=1; end	% step
%if (nargin<4), ws=5; end	% winsize
 
[n,m]=size(sq);


% [patt, npatt, scate] = sitepattern(sq);
% nop=length(npatt);
% lpatt=zeros(1,nop);
% fprintf (['total sites=%d\n'], m)
% fprintf (['total patterns=%d\n'], nop)
% for (k=1:nop),
%     fprintf (['calculating likelihood for pattern %d/%d ...\n'], k,nop)
%     lpatt(1,k)=treelike(patt(:,k),tr,md);
% end
% siteL=lpatt(scate);

[totL,siteL] = treelike(tr,sq,md,1,0);
totL=sum(siteL);

disp(sprintf('Tree -ln Likelihood (total) = %f',totL))

minws=5;
halplen=round(m/2);
Q=[];
%Q=zeros(1,halplen-minws+1);

%sp=1;
%slidingfun(@i_platoQ,siteL,minws,sp,totL,restlen)
%pause

for (ws=minws:halplen),
	restlen=length(siteL)-ws;
	sp=1;
	q=slidingfun(@i_platoQ,siteL,ws,sp,totL,restlen);
	Q=[Q;q];
	%pause
	%Q(ws-4)=max(find(q==max(q)));
end

[x,ix]=max(Q);



%v=slidingfun(@treelike,sq,ws,st,tr,md);
%if (nargout<1),
%	plot([1:length(q)],q)
%end



function [q] = i_platoQ(winL,totL,restlen)
%This value is divided by a
%similar measure of the average log likelihood of the rest
%of the sequence (of length n) to give the final measure,
%Q. This gives a measure of the average likelihood of the
%window with respect to the rest of the sequence
%This gives a measure of the average likelihood of the
%window with respect to the rest of the sequence that is
%consistent over different window sizes.

avgRest=(totL-sum(winL))/restlen;
q=mean(winL)/avgRest;




%		for(size=minWinSize;size<half_plus_one;size++){
%			for(sp=0;sp<(numBases+1-size);sp++){
%				temp=0.0;			
%				for(j=sp;j<(sp+size);j++)
%					temp+=simLik[j];
%				temp=(temp/size)/((likelihood-temp)/(numBases-size));
%				if(temp>winLikList[size-minWinSize][i])
%					winLikList[size-minWinSize][i]=temp;
%			}
%		}
