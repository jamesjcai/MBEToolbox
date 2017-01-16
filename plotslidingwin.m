function plotslidingwin(aln,analysistype,winsize)
%PLOTSLIDINGWIN - Performs sliding window analysis on a nucleotide sequence
%
% Syntax: plotslidingwin(aln,'analysistype',winsize)
%
% Inputs:
%    aln            - Distance matrix
%    analysistype   - 'GC'|'AG'|'GCDeviation'|'ATDeviation'
%    winsize        - Window size
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



warning off MATLAB:divideByZero;

Seq=aln.seq(1,:);
disp('Working on the first sequence in alignment.');
if nargin < 2 | isempty(analysistype), analysistype = 'GC'; end
if nargin < 3 | isempty(winsize) |  isnan(winsize), winsize = 0; end


switch (analysistype)
    case ('GC')
	if (winsize==0), winsize=120; end
	GC=Seq==2|Seq==3;

	plot(slidingavg(GC,winsize))
	info = ['GC Content (%) Window size ', num2str(winsize)];
	title(info);
	axis([1 length(GC) min(GC)*1.1 max(GC)*1.1]);
	xlabel('Base (bp)'); ylabel('GC (%)');
    case ('AG')
	if (winsize==0), winsize=120; end
	AG=Seq==1|Seq==3;

	plot(slidingavg(AG,winsize))
	info = ['AG Content (%) Window size ', num2str(winsize)];
	title(info);
	axis([1 length(AG) min(AG)*1.1 max(AG)*1.1]);
	xlabel('Base (bp)'); ylabel('AG (%)');
    case ('GCDeviation')
	if (winsize==0), winsize=30; end
	C=Seq==2; G=Seq==3;
	SC = slidingavg(C,winsize);
	SG = slidingavg(G,winsize);
	GCDev = (SG-SC)./(SG+SC);

	plot(GCDev)
	info = ['GC Deviation (G-C)/(G+C) Window size ', num2str(winsize)];
	title(info);
	axis([1 length(GCDev) min(GCDev)*1.1 max(GCDev)*1.1]);
	xlabel('Base (bp)'); ylabel('GC Deviation');

    case ('ATDeviation')
	if (winsize==0), winsize=30; end
	A=Seq==1; T=Seq==4;
	SA = slidingavg(A,winsize);
	ST = slidingavg(T,winsize);
	ATDev = (SA-ST)./(SA+ST);

	plot(ATDev)
	info = ['AT Deviation (A-T)/(A+T) Window size ', num2str(winsize)];
	title(info);
	axis([1 length(ATDev) min(ATDev)*1.1 max(ATDev)*1.1]);
	xlabel('Base (bp)'); ylabel('AT Deviation');
end