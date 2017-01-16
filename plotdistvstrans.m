function plotdistvstrans(aln)
%PLOTDISTVSTRANS - Plot genetic distances (JC69) vs. transitions and transversions
%
% Syntax: plotdistvstrans(aln)
%
% Inputs:
%    aln      - Alignment structure
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




[n,m] = size(aln.seq);
[P,Q]=countseqpq(aln);
P = P./m;	% transition
Q = Q./m;	% transversion
%D = dn_jc(aln);
D = dn_hky(aln);

l=1;
for (i=1:n)
    for (j=i:n)
	if (i ~= j)
		dp(l,1)=D(i,j);
		dp(l,2)=P(i,j);
		dq(l,1)=D(i,j);
		dq(l,2)=Q(i,j);
	        l=l+1;
	end
    end
end


% graphical output

plot(dp(:,1), dp(:,2),'g+','MarkerSize',5);
hold on;
plot(dq(:,1), dq(:,2),'v','MarkerSize',5);
info = 'Transitions & Transversions vs. Distance';
title(info);
xlabel('Distance (HKY)'); ylabel('Transitions & Transversions');
legend('Transition','Transversion');
hold off;