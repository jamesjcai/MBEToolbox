function [rateraw,ratenorm]=rateofsites_ml(seq,tree,model,noise)
%RATEOFSITES_ML - Maximum-likelihood estimation of site-specific evol. rate
%
% Syntax: [rateraw,ratenorm]=rateofsites_ml(seq,tree,model)
%
% Inputs:
%    seq        - Sequences
%    tree       - Tree, a string like '((1:0.23,2:0.23):0.45,3:0.10);')
%    model      - model, a struct created by MODEL* functions
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
% See also: RATEOFSITES_EB

% Molecular Biology & Evolution Toolbox, (C) 2006-2010
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://www.hku.hk/jamescai/mbetoolbox
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<4), noise=0; end
if (nargin<3), error('Need at least 3 parameters.'); end
if (isstruct(seq)), seq=seq.seq; end

m=size(seq,2);
rateraw=zeros(1,m);
[patt, npatt, scate] = sitepattern(seq);
nop=length(npatt);
rpatt=zeros(1,nop);
fprintf (['total sites=%d\n'], m)
fprintf (['total patterns=%d\n'], nop)
for (k=1:nop),
    fprintf (['optimizing rate for pattern %d of %d ...\n'], k, nop)
    rpatt(1,k)=optimsiterate(patt(:,k),tree,model);
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
	    fprintf (['\n#pos\tseq\tscore\tnormalised score\n']);
	for (k=1:m)
	    fprintf (['%d\t%s\t%5.4f\t%10.4f\n'],k, X(seq(1,k)), rateraw(1,k), ratenorm(1,k));
	end
	disp('#Average of normalised score = 0.0000')
	disp('#Standard deviation of normalised score = 1.0000')
end