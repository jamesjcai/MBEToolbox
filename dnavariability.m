function [res]=dnavariability(s,showit)
%DNAVARIABILITY - DNA variability
%This function estimates variability in ech column of the alignment as an
%entropy function of the nucleotide variantion obseved in this column using the
%following equation:
%
%$$Var = \sum\limits_{i = A,C,G,T}^{} {\frac{{n_i}}{N}} \ln \frac{{n_i}}{N}$$
%
% Reference:
%  Proutski V, Holmes E.
%  SWAN: sliding window analysis of nucleotide sequence variability.
%  Bioinformatics. 1998 Jun;14(5):467-8.
%  PMID: 9682061
%
% Syntax: [res]=dnavariability(s,showit)
%
% Inputs:
%    s            - Sequences
%    showit       - If showit=1, then plot window analysis
%
% Outputs:
%    res      - DNA variability
%
% Example:
%
% >>s=[4 4 1 1;4 4 1 1;4 4 1 1;4 2 1 1;4 2 1 1;2 1 1 1;2 1 1 1;2 1 1 3;2 3 1 3;2 3 3 3];
% >>[res]=dnavariability(s,showit)
%
% See also: COMPCOMP

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $




if (nargin<2),
    showit=0;
end
[n,m]=size(s);
res=zeros(1,m);
for (k=1:m),
     col=s(:,k);
     res(1,k)=i_colvariability(col,n);
 end


 res=-1*res;
if (showit==1),
 	winsize=30;
	dnav = slidingavg(res,winsize);
    resavg=mean(dnav);
    resstd=std(dnav);
	plot(dnav);
	info = ['DNA Variability - window size ', num2str(winsize)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Base (bp)'); ylabel('DNA Variability');

    hold on
    plot([1:m],ones(1,m)*resavg,'r')
    plot([1:m],ones(1,m)*(resavg+resstd),'g--')
    plot([1:m],ones(1,m)*(resavg-resstd),'g--')
	legend('', 'mean','std') ;
    hold off
end



function [res]=i_colvariability(col,n)
res=0;
for (k=1:4)
     x=sum(col==k);
     if (x>0),
         y=x./n;
         res=res+y*log(y);
     end
 end

