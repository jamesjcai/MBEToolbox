function [lnL] = siteratelike(s,tr,m,rate)
%SITERATELIKE - Negative gamma log-likelihood function.

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


[Q] = composeQ(m.R, m.freq);
[W,L] = eig(Q); V = inv(W);
freq = m.freq;

%rate = 1;
nalpha = length(m.R);

[tree, n_otu, t, nm] = parsetree(tr);
[n_otu2, n_site]=size(s);
if ~(n_otu==n_otu2)
	error('Inconsistent tree and alignment. Need an rooted tree.');
else
	inlineL=mbelfcreator(tree,n_otu);
end

% t = ones(1,n_brc);
n_brc = 2*n_otu - 3;
if ~(length(t)-2==n_brc),
	error('Inconsistent tree and alignment. Need an rooted tree.');
else
	P = zeros(n_brc*nalpha,nalpha);
	for j=1:n_brc
	    P((j-1)*nalpha+1:j*nalpha,:)=W*expm(L*t(j)*rate)*V; 
	end
end

lnL=0;
for (k=1:n_site),
	site = s(:,k);
	Lhi = inline(inlineL,'P','n_otu','freq','s','nalpha');
	lnL = lnL+log(Lhi(P,n_otu,freq,site,nalpha));
end
