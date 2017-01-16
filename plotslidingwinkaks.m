function [dS,dN,dS_cum,dN_cum] = plotslidingwinkaks(aln,winsize,option)
%PLOTSLIDINGWINKAKS - Plots sliding Ka and Ks curves
%
% Syntax: [dS,dN,dS_cum,dN_cum] = plotslidingwinkaks(aln,winsize,option)
%
% Inputs:
%    aln       - Alignment structure
%    winsize   - Window size
%    option    - 1 - plot sliding window of Ka and Ks; 2 - plot cumulative Ka and
%                Ks curves; 3 - plot enhanced sliding window of Ka and Ks
%
% Outputs:
%    dS       - Vector containing synonymous substitution rate (Ks) along sequences
%    dN       - Vector containing nonsynonymous substitution rate (Ka) along sequences
%    dS_cum   - Cumulative dS
%    dN_cum   - Cumulative dN
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


if nargin < 2 | isempty(winsize) | isnan(winsize), winsize = 20; end
if nargin < 3, option = 1; end

if ~(isvalidaln(aln,'CODING')),
	error ('ERROR: Not coding seq')
end
if (hasgap(aln))
	disp ('Gaps in sequences have been removed.')
	aln=rmcodongaps(aln);
end
[n,m3]=size(aln.seq);

if (n<2)
    error ('ERROR: at least two sequences required');
elseif (n>2)
    disp ('WARNING: only first two sequences were used.')
    aln.seq=aln.seq([1 2],:);
end


m=m3/3;
[CodonSeq]=codonise64(aln.seq);
[S,N]=getsynnonsynsites(aln.geneticcode);
[ns,na]=getsynnonsyndiff(aln.geneticcode);
vector_s = sum(S(CodonSeq))./2;
vector_n = sum(N(CodonSeq))./2;
vector_nst = zeros(1,m);
vector_nat = zeros(1,m);
for (k=1:m),
	vector_nst(1,k) = ns(CodonSeq(1,k),CodonSeq(2,k));
	vector_nat(1,k) = na(CodonSeq(1,k),CodonSeq(2,k));
end
vec_s = slidingsum(vector_s,winsize);
vec_n = slidingsum(vector_n,winsize);
vec_nst = slidingsum(vector_nst,winsize);
vec_nat = slidingsum(vector_nat,winsize);
[dS,dN] = NG86_Correction(vec_s,vec_n,vec_nst,vec_nat);

if (length(dN)~=m|length(dS)~=m)
      error('Something wrong!')
end


dS_cum = cumsum(dS);
dN_cum = cumsum(dN);


switch (option)
    case (1)
    	plot([1:m],dS,[1:m],dN);
    	if (max(max(dN),max(dS))>=2.5)
    		axis([1 m 0 2.5]);
    	else
    		axis([1 m 0 max(max(dN),max(dS))*1.1]);
    	end
       	ylabel('Substitution rate');
 	    legend('syn','nonsyn');
	    xlabel('Codon site');
    case (2)
    	plot([1:m],dS_cum,[1:m],dN_cum);
    	axis([1 m 0 max(dS_cum)*1.1]);
    	ylabel('Cumulative substitution number per site');
	    xlabel('Codon site');
    case (3)
    	dS_amp = i_zamplify(dS_cum);
    	dN_amp = i_zamplify(dN_cum);

        % USELESS NORMALISATION
        % dS_amp=(dS_amp-mean(dS_amp))./std(dS_amp);
        % dN_amp=(dN_amp-mean(dN_amp))./std(dN_amp);
        % dS_amp=diff(dS_amp); dN_amp=diff(dN_amp);
    	% plot([1:length(dS_amp)],dS_amp,[1:length(dS_amp)],dN_amp);

    	plot([1:m],dS_amp,[1:m],dN_amp);
       hold on
        highlightpeaks(dN_amp);
        highlightpeaks(dS_amp);
       hold off

       	axis([1 m min([dS_amp dN_amp])*1.1 max([dS_amp dN_amp])*1.1]);
    	ylabel('Z'' transformed substitution rate');
	    xlabel('Codon site');
    case (5)
    	plot([1:m],dN);
    end

	%info = ['Method: NG86; Window size: ', num2str(winsize)];
	%title(info);
	% legend('synonymous','nonsynonymous') ;



function highlightpeaks(X)
        Xmax = find(diff(diff(X)<0)==1)+1;
        Xmin = find(diff(diff(X)>0)==1)+1;
        for (k=1:length(Xmax))
            x=Xmax(k); y=X(x);
            plot(x,y,'r.')
        end
        for (k=1:length(Xmin))
            x=Xmin(k); y=X(x);
            plot(x,y,'r.');
        end



function [dS,dN] = NG86_Correction(s,n,nst,nat)

pS=nst./s; pN=nat./n;
dS = 1-(4/3)*pS;
dN = 1-(4/3)*pN;

dS(find(dS~=1&dS>=0))=-0.75*log(max(eps,dS(find(dS~=1&dS>=0))));
dS(find(dS==1))=0;
dS(find(dS<0))=-1;

dN(find(dN~=1&dN>=0))=-0.75*log(max(eps,dN(find(dN~=1&dN>=0))));
dN(find(dN==1))=0;
dN(find(dN<0))=-1;




function [Z1] = i_zamplify(Z)
[n,m]=size(Z);
Xaxis = 1:m;
Yaxis = Z;
[b,m,E]=lsquare(Xaxis',Yaxis');	% m is slope
Z1=Z-m*Xaxis;