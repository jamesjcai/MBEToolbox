function [Q] = plato(sq,tr,md)
%PLATO - Detecting recombination with window methods
%Partial Likelihoods Assessed Through Optimisation
%
% Syntax: plato(sq,tr,md)
%
% Inputs:
%    sq     - Seq
%    tr     - Tree string
%    md     - Model
%
% Outputs:
%    Q    - Raw Q value by PLATO method
%
% REF: Grassly, Holmes (Mol. Biol. Evol. 14)
%
% See also:

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $



if (isstruct(sq)), sq=sq.seq; end
[n,m]=size(sq);
[totL,siteL] = treelike(tr,sq,md,1,0);
totL=sum(siteL);

disp(sprintf('Tree -ln Likelihood (total) = %f',totL))

minws=5;
halplen=round(m/2);
Q=zeros(halplen-minws+1,m-1);

for (ws=minws:halplen),      % window size from 5 to half length of seq
	restlen=m-ws;
	q=slidingfun(@i_platoQ,siteL,ws,1,totL,restlen);
	Q(ws-4,:)=q;
end


if (nargout<1),
      plot(Q');
      axis([1 m min(Q(:))*1.1 max(Q(:))*1.1]);
      title('PLATO - Detecting recombination with window methods');
      xlabel('Base (bp)'); ylabel('Raw Q(i,j)');
      disp('Note: the Q(i,j) values reported have not been normalized.')
end



function [q] = i_platoQ(winL,totL,restlen)
%This q value gives a measure of the average likelihood of the window with respect 
%to the rest of the sequence that is consistent over different window sizes.

avgRest=(totL-sum(winL))/restlen;
q=mean(winL)/avgRest;
