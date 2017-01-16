function [aln2]=rmcodongaps(aln)
%RMCODONGAPS - Remove codons with gaps
%
% Syntax:  [aln2]=rmcodongaps(aln)
%
% Inputs:
%    aln    - Alignment structure
%
% Outputs:
%    aln2   - Alignment structure
%
%
% See also: HASGAP

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (isstruct(aln)),
	if ~(isvalidaln(aln,'CODING')),
		error ('ERROR: Not coding seq')
	end


	[aln2] = copyalnheader(aln);
	if ~(hasgap(aln))
		aln2.seq=aln.seq;
		return;
	end

	S=aln.seq;
	[n,m] = size(S);
	S2=[];
	for (k=1:3:m),
		X=S(:,(k:k+2));
		if ~(sum(sum(X>4))>0)
			S2=cat(2,S2,X);
		end
	end
	aln2.seq=S2;

else
	seq=aln;
	S=seq;
	[n,m] = size(S);
	S2=[];
	for (k=1:3:m),
		X=S(:,(k:k+2));
		if ~(sum(sum(X>4))>0)
			S2=cat(2,S2,X);
		end
	end
	aln2=S2;
end