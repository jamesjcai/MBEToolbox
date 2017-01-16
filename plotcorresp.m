function plotcorresp(aln,option)
%PLOTCORRESP - Performs correspondence analysis and plot result
%
% Syntax: plotcorresp(aln,'option')
%
% Inputs:
%    aln      - Alignment structure
%    option   - (optional) 'rscu'|'codon', COA on RSCU or codon usage
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


if(nargin<2),
   option='rscu';
end

switch (option)
    case ('codon')
	[K] = codonusage(aln);
    case ('rscu')
	[K] = rscu(aln,'notall');
   otherwise
	[K] = rscu(aln,'notall');
end

[m,n] = size(K);
r = min(m-1,n);     % max possible rank of K
avg = mean(K);		%  Computes and prints means
centerx = (K - avg(ones(m,1),:));	%  Computes deviations from Means

	 % centerx./sqrt(m-1)		% Computes variance-covariance
[U,latent,pc] = svd(centerx./sqrt(m-1),0);	% Singular value decomposition
						% of varcov for
						% eigenvalues L and
						% eigenvectors E
% pc = eigenvectors E
score = centerx*pc;				% PC scores

plot(score(:,1),score(:,2),'+b')  % : selects all rows
                                    % and plots
                                    % column 1 K 2
                                    % of PC Scores
xlabel('FIRST AXIS')
ylabel('SECOND AXIS')
title('Correspondence Analysis of Codon Usage')
% title('Plot of PC2 against PC1')

num=int2str([1:m]');
text(score(:,1),score(:,2),num)
grid
hold off

for (k=1:m),
	fprintf(['%3d = %s\n'], k, aln.seqnames{k});
end
