function plotkavsks(aln)
%PLOTKAVSKS - Plot log(Ka) vs. log(Ks)
%
% Syntax: plotkavsks(aln)
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



aln2=rmgaps(aln);
[n,m] = size(aln2.seq);


ButtonName=questdlg('What algorithm you want to use?', ...
                    'Select algorithm', ...
                    'NG86','Li85/93','ML','NG86');
switch upper(ButtonName),
    case 'NG86',
     [dS,dN,dN_dS,St,Nt,nst,nat,VdS,VdN]=dc_ng86(aln2);
     case 'LI85/93',
     [dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln2);
     case 'ML',
     [dS,dN,dN_dS]=dc_ml(aln2);
     otherwise
     [dS,dN,dN_dS,St,Nt,nst,nat,VdS,VdN]=dc_ng86(aln2);
end






dr=zeros(0.5*(n^2-n),2);

ct=0;
for (i=1:n)
    for (j=i:n)
	if (i ~= j)
		if ~(isnan(dS(i,j))|isnan(dN(i,j)))
			ct=ct+1;
			dr(ct,1)=dS(i,j);
			dr(ct,2)=dN(i,j);
		end
	end
    end
end

da=dr(:,1);
db=dr(:,2);
%loglog(da, db, 'o')
plot(da,db,'o','MarkerSize',5);
infostr=sprintf('Ka vs. Ks, (%d seqences, %d points)',n,0.5*(n^2-n));
title(infostr);
xlabel('Ks'); ylabel('Ka');
axis equal
%axis([0 2 0 2]);
% grid;

hold on
mx=max(max(da),max(db));
mi=min(min(da),min(db));
plot([mi mx],[mi mx])
axis([mi mx mi mx]);
hold off
