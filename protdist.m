function [D]=protdist(aln,mdlname)

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2)
	mdlname='jtt';
end

[n,m] = size(aln.seq);
D = zeros(n);

switch (lower(mdlname))
    case {'jtt'}
	 model=modeljtt;
    case {'dayhoff'}
	model=modeldayhoff;
    case {'wag'}
	model=modelwag;
    case {'lg'}
	model=modellg;
   otherwise
	model=modeljtt;
end

%	fprintf(['\n']);
for i=1:n
	fprintf(['%s '],aln.seqnames{i});
for j=i:n
if (i~=j)
	s1=aln.seq(i,:);
	s2=aln.seq(j,:);
	D(i,j) = optimseqpairlikeli(model,s1,s2);
	fprintf(['.']);
end
D(j,i) = D(i,j);
end
	fprintf(['\n']);
end