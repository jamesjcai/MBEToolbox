function [p] = compequtest(s,freq)
%COMPEQUTEST - Compositional equilibrium test
%
%Maximum likelihood phylogenetic inference is based on the assumption that the set of
%homologous sequences is in compositional equilibrium. This is postulating the existence of
%selective forces that keep the nucleotide frequencies in a sequence constant over time. As the
%violation of assumptions may produce the wrong results, one should always check their
%validity in advance of a phylogenetic analysis. --- <<TREEFinder Mannual>>

%For each sequence in the chosen alignment, the utility calculates the p-value that its nucleotide
%distribution does not differ from the expectation. Their p-values should be all well above 0.05.
%Otherwise, it becomes improbable that the sequence has the expected base composition and
%one should not assume equilibrium.
%
% REF: TREEFINDER MANUAL - Version of November 2004 - Gangolf Jobb

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if isstruct(s)
      names=s.seqnames;
      s=s.seq;
else
      names={};
end


if (nargin<2)
	freq=estimatefreq(s);
end
disp('==============================')
disp('COMPOSITIONAL EQUILIBRIUM TEST')
disp('==============================')
disp(' ')
disp('This is a chi-square test comparing the nucleotide composition of');
disp('each sequence to that of the whole data set.');
disp(' ');
disp('Expectation:');
disp(' ');
fprintf(['  A: %2.5f\n'],freq(1));
fprintf(['  C: %2.5f\n'],freq(2));
fprintf(['  G: %2.5f\n'],freq(3));
fprintf(['  T: %2.5f\n'],freq(4));
disp(' ');
disp('Nucleotide counts A,C,G,T: ');

[n,m]=size(s);
M=[sum(s==1,2),sum(s==2,2),sum(s==3,2),sum(s==4,2)];
p=zeros(1,n);
ExpectM = m.*freq;


disp(' ');
for (k=1:n),
      fprintf(['%8s :%5d,%5d,%5d,%5d\n'], i_names(k,names), M(k,:));
end
disp('    -----------------------------');
fprintf(['%8s :%5d,%5d,%5d,%5d\n'], ['Total'], sum(M,1));
disp(' ');
disp('P-Values:');
disp(' ');
disp('WARNING: This test is statistically unreliable!');

for (k=1:n),
	x=sum(((M(k,:)-ExpectM).^2)./ExpectM);
	if ~(x>1400) % otherwise p=0.0
	      p(1,k)=chi2calc(x,4-1);
	end
	if (p(1,k)>0.05)
		fprintf(['%8s : %2.5f\n'], i_names(k,names), p(1,k));
	else
		% give warning of possible unequilibrium composition
		fprintf(['%8s : %2.5f <==!!!\n'], i_names(k,names), p(1,k));
	end
end
disp(' ');
disp('All p-values should be well above 0.05 at compositional equilibrium.');


function [txt] = i_names(k,names)
if (isempty(names))
	txt=['seq',num2str(k)];
else
	txt=char(names(k));
end

