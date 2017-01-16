function [aln2]=rmgaps(aln)
%RMGAPS - Remove gaps in alignment
%
% Syntax:  [aln2]=rmgaps(aln)
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



[n,m] = size(aln.seq);
[aln2] = aln;

codes = zeros(1,m);

switch (aln.seqtype)
    case (3)

    	for i = 1:m
		for j = 1:n
		if ((aln.seq(j,i) == '-') | (aln.seq(j,i) == ' ') |(aln.seq(j,i) >= 21))
		codes(i) = 1;
				break;
			end
		end
	end

    otherwise

    	for i = 1:m
		for j = 1:n
		if ((aln.seq(j,i) == '-') | (aln.seq(j,i) == ' ') |(aln.seq(j,i) == 5))
		codes(i) = 1;
				break;
			end
		end
	end

end



for i = 1:m
	j = m - i + 1;
   if codes(j) == 1
	aln2.seq(:,j) = [];
   end
end







% method=2;
% switch (method)
%     case (1)
% 	[n,m] = size(aln.seq);
% 	picker=aln.seq(n,:)~=5;
% 	[aln2] = copyalnheader(aln);
% 	for i=1:n
% 		aln2.seq(i,:) = aln.seq(i,picker);
% 	end
%   case (2)
% 	[n,m] = size(aln.seq);
%	[aln2] = aln;
%	codes = zeros(1,m);
%	for i = 1:m
%		for j = 1:n
%		if ( (aln.seq(j,i) == '-') | (aln.seq(j,i) == ' ') |(aln.seq(j,i) == 5))
%		codes(i) = 1;
%				break;
%			end
%		end
%	end
%
%
%	for i = 1:m
%		j = m - i + 1;
%	   if codes(j) == 1
%		aln2.seq(:,j) = [];
%	   end
%	end
%	end