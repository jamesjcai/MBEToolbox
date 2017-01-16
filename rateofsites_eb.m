function [rateraw,ratenorm]=rateofsites_eb(seq,tree,model,ncate,noise)
%RATEOFSITES_EB - Empirical Bayesian estimation of site-specific evol. rate
%
% Syntax: [rateraw,ratenorm]=rateofsites_eb(seq,tree,model,ncate)
%
% Inputs:
%    seq        - Sequences
%    tree       - Tree, a string like '((1:0.23,2:0.23):0.45,3:0.10);')
%    model      - Model, a struct created by MODEL* functions
%    ncate      - Number of categories in a gamma distribution (default=16)
%
% Outputs:
%    rateraw    - Evolutionary rates of sites
%    ratenorm   - Normalised evolutionary rates of sites
%
% Reference:
% Mayrose I, Graur D, Ben-Tal N, Pupko T. Comparison of site-specific
% rate-inference methods for protein sequences: empirical Bayesian methods are
% superior. Mol Biol Evol. 2004 Sep;21(9):1781-91
%
% See also: RATEOFSITES_ML

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://www.hku.hk/jamescai/mbetoolbox
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<5), noise=0; end
if (nargin<4), ncate=16; end
if (nargin<3), error('Need at least 3 parameters.'); end

if (isstruct(seq)), seq=seq.seq; end

disp('Estimating global rate...')
[gratev,gratevn]=rateofsites_ml(seq,tree,model);  % global rate
grate=mean(gratev);

fprintf('global rate, grate=%2.4f\n',grate)
disp('Estimating alpha of gamma distribution...')
x=gamfit(gratev);
alpha=x(1);
gd=gammadistrib(ncate,alpha);
fprintf('alpha=%2.4f, ncate=%d (default=16)\n',alpha,ncate)


[n,m]=size(seq);
rateraw=zeros(1,m);
[patt, npatt, scate] = sitepattern(seq);
nop=length(npatt);
rpatt=zeros(1,nop);
fprintf (['total sites=%d\n'], m)
fprintf (['total patterns=%d\n'], nop)

for (k=1:nop),
    % fprintf (['estimating rate for pattern %d ...\n'], k)
    rpatt(1,k)=i_rateofsite_eb(patt(:,k),tree,model,gd,grate);
end

rateraw=rpatt(scate);
if (nargin>1),
	ratenorm=(rateraw-mean(rateraw))./std(rateraw);
end

if ((nargout<1)||(noise)),
	[NT,AA] = seqcode;
	if (length(model.R)==4),
	    X=NT;
	elseif(length(model.R)==20),
	    X=AA;
	end
	    fprintf (['\n#Pos\tSeq\tRaw score\tNml. score\n']);
	for (k=1:m)
	    fprintf (['%d\t%s\t%5.4f\t%10.4f\n'],k, X(seq(1,k)), rateraw(1,k), ratenorm(1,k));
	end
	disp('#Average of normalised score = 0.0000')
	disp('#Standard deviation of normalised score = 1.0000')
end



function [er]=i_rateofsite_eb(seq,tree,model,gd,grate)
	%er - Expected r
	ri=gd.rate;
	pri=gd.prob;    % P(r_i) is the prior distribution on the rate = 1/ncate
	ncate=length(ri);
	rcate=zeros(1,ncate);

	% issite=1;
	for (k=1:ncate),
	%    rcate(1,k) = -1*siteratelike(seq,tree,model,grate*ri(k)*pri(k));
	    xlnL=treelike(seq,tree,model,grate*ri(k)*pri(k));
	    rcate(1,k) = exp(xlnL);  % siteratelike
	    % fprintf ('fhK=%5.5f\n',rcate(1,k))
	end

	%Posterior probabilities for site classes
	xx=max(sum((rcate.*pri)),eps);
	pr=(rcate.*pri)./xx;
%	er=sum(pr.*rcate);   9/29/2010 %Thank Jorge Enrique Moreira Broche <jmoreira@uclv.edu.cu> for correcting this error.
	er = sum(pr.*(grate.*ri));
	%fprintf ('fhK=%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5d\n',log(xx),xx,128*xx,er,find(rcate==max(rcate)))
