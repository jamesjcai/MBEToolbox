function [v] = slidingfun(funfcn,aln,winsize,step,varargin)
%SLIDINGFUN Scalar bounded nonlinear function minimization.  
%   USAGE: [v]=slidingfun(funfcn,aln,winsize,step,varargin)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

funfcn = fcnchk(funfcn,length(varargin));
%if ~isa(funfcn,'function_handle'), ...
if (nargin<3), winsize=20; end
if (nargin<4), step=1; end

[n,m]=size(seq);
if (winsize>m|step>m), error('Please select appropriate winsize and step'); end
if (winsize<step), error('Please select appropriate winsize and step'); end
if (winsize<1|step<1), error('Please select appropriate winsize and step'); end

sp = 1:step:m;
ep = sp + (winsize-1);
ep(find(ep>m))=m;
nwin = length(sp)-1;
v=zeros(1,nwin);


%h = waitbar(0,'Please wait...');
for (k=1:nwin),
      v(k)=feval(funfcn,seq(:,sp(k):ep(k)),varargin{:});
%      waitbar(k/nwin)
end
%close(h) 
if (nargout<1),
	plot(v);
	res=v;
	resavg=mean(res);
	resstd=std(res);
	score=sum(res>=resavg+resstd)+sum(res<=resavg-resstd);
	[n,m]=size(res);
	info = ['Data - window size ', num2str(winsize), ' step ', num2str(step)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Site Position'); ylabel('Statistic');
	    hold on
	    plot([1:m],ones(1,m)*resavg,'r')
	    plot([1:m],ones(1,m)*(resavg+resstd),'g--')
	    plot([1:m],ones(1,m)*(resavg-resstd),'g--')
	legend('statistic', 'mean','std') ;
    hold off
end




function [v] = i_old_slidingfun(funfcn,aln,winsize,step,varargin)
%SLIDINGFUN Scalar bounded nonlinear function minimization.  
%   USAGE: [v]=slidingfun(funfcn,aln,winsize,step,varargin)

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

funfcn = fcnchk(funfcn,length(varargin));
%if ~isa(funfcn,'function_handle'), ...
if (nargin<3), winsize=20; end
if (nargin<4), step=1; end

[n,m]=size(seq);
if (winsize>m|step>m), error('Please select appropriate winsize and step'); end
if (winsize<step), error('Please select appropriate winsize and step'); end
if (winsize<1|step<1), error('Please select appropriate winsize and step'); end

sp = 1:step:m;
ep = sp + (winsize-1);
ep(find(ep>m))=m;
nwin = length(sp)-1;
v=zeros(1,nwin);


%h = waitbar(0,'Please wait...');
for (k=1:nwin),
      v(k)=feval(funfcn,seq(:,sp(k):ep(k)),varargin{:});
%      waitbar(k/nwin)
end
%close(h) 
showit=1;
if (showit==1),
	plot(v);
res=v;
resavg=mean(res);
resstd=std(res);
score=sum(res>=resavg+resstd)+sum(res<=resavg-resstd);
[n,m]=size(res);

	info = ['Data - window size ', num2str(winsize)];
	title(info);
	%axis([1 length(dnav) min(dnav)*1.1 max(dnav)*1.1]);
	xlabel('Nucleotide position'); ylabel('D');
    hold on
    plot([1:m],ones(1,m)*resavg,'r')
    plot([1:m],ones(1,m)*(resavg+resstd),'g--')
    plot([1:m],ones(1,m)*(resavg-resstd),'g--')
	legend('D', 'mean','std') ;
    hold off
end
