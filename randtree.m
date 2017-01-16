function rt = randtree(n)
%RANDTRE - Generates a random tree with n OTUs
%
% Adopted from PHYLLAB toolbox v1.1

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<1), n=3; end
if (n<3), error('Tree should contain at least 3 elements.'); end


otu = randperm(n);
% unite the first pair (in any tree is at least one pare!!)
brlen=rand;
rt = sprintf('(Seq%d:%f,Seq%d:%f)', otu(1), rand, otu(2), brlen);

% process the others
flag = 0;
for i=3:n
	if flag
		flag = 0;
	else

		if rand > 0.5
			% it's a standalone otu
			brlen=rand;
			rt = sprintf('(%s:%f,Seq%d:%f)',rt,rand,otu(i),brlen);
		else

			% it's a paired otu
			if i+1 <= n
				brlen=rand;
				rt = sprintf('(%s:%f,(Seq%d:%f,Seq%d:%f):%f)', rt, rand, otu(i), rand, otu(i+1), brlen, rand);
				flag = 1;
			else
				brlen=rand;
				% it can't be paired, we running out of availible otu's
				rt = sprintf('(%s:%f,Seq%d:%f)',rt,rand,otu(i),brlen);
			end;
		end;
	end;
end;
