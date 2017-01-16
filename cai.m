function [cai] = cai(aln,refindex)
%CAI - Calculates the Codon Adaptation Index (cai)
%Codon adaptation index (cai) (Sharp PM & Li WH 1986) has previously been
%described for E.coli and yeast. Codon usage in a reference set of highly
%expressed genes is used to estimate the "relative adaptiveness", w, of each
%codon. w is calculated as the frequency of use of that codon (in the reference
%set) relative to the frequency of the optimal codon for that amino acid. The
%cai for any gene is estimated as the geometric mean of the w values
%corresponding to each of the codons in that gene (excluding Met, Trp and
%termination codons.) The maximum possible value for cai is 1.0, when only
%optimal codons are used.
%
% Syntax: [cai] = cai(aln,refindex)
%
% Inputs:
%    aln      - Alignment structure
%    refindex - (optional) Index of reference codon usage set (default=1)
%               1 - Escherichia coli (Sharp PM & Li WH 1986)
%               2 - Bacillus subtilis
%               3 - Saccharomyces cerevisiae (Sharp and Cowe 1991 Yeast 7:657-678)
%
% Outputs:
%    cai     - Codon Adaptation Index value
%
% Examples:
%
% >> [cai] = cai(aln);
% >> [cai] = cai(aln,1);
% % Calculate cai using B. subtilis reference codon usage set
% >> [cai] = cai(aln,2);
%
% See also: RSCU

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if ~isvalidaln(aln,'CODING'), error('Not coding seq'); end
if hasgap(aln)
	disp ('WARNING: Sequences cannot contain gaps. Gaps will be removed by using REMOVEGAPS.')
	aln=rmgaps(aln);
end
if nargin<2, refindex=1; end


% To calculate cai a reference set of highly expressed genes
% must be selected:
% 1 - Escherichia coli, Sharp PM & Li WH (1986)
% 2 - Bacillus subtilis
% 3 - Saccharomyces cerevisiae, Sharp and Cowe (1991) Yeast 7:657-678
W=[1.000 1.000 0.135;
1.000 1.000 1.000;
0.253 0.097 1.000;
0.051 0.417 0.053;
0.076 0.867 0.012;
1.000 0.033 1.000;
0.099 0.200 0.006;
0.965 1.000 0.921;
0.004 0.435 1.000;
0.410 0.208 0.031;
0.002 0.022 0.003;
0.085 0.125 0.021;
0.003 0.071 0.003;
1.000 1.000 1.000;
1.000 1.000 1.000;
0.185 0.500 0.823;
0.124 1.000 1.000;
1.000 0.083 1.000;
1.000 0.214 0.007;
0.291 1.000 0.245;
0.135 0.714 1.000;
0.012 0.071 0.009;
1.000 0.143 0.002;
0.070 1.000 0.047;
0.004 0.022 0.002;
0.356 0.609 0.002;
0.004 0.043 0.002;
1.000 1.000 0.137;
0.007 0.500 0.039;
0.037 0.143 0.003;
1.000 0.071 0.003;
0.042 0.857 0.006;
1.000 1.000 1.000;
1.000 1.000 1.000;
0.259 0.412 0.016;
0.434 0.417 0.554;
0.586 0.275 0.015;
0.122 0.025 0.316;
0.424 0.125 0.001;
1.000 1.000 1.000;
0.010 1.000 0.002;
0.724 0.773 0.020;
0.019 0.045 0.004;
1.000 0.955 1.000;
0.495 0.750 0.002;
0.066 0.188 0.831;
0.221 0.438 0.018;
1.000 1.000 1.000;
0.000 0.000 0.000;
1.000 1.000 1.000;
0.000 0.000 0.000;
0.239 0.500 0.071;
0.077 0.458 0.036;
0.744 0.021 0.693;
0.017 0.021 0.005;
1.000 1.000 1.000;
0.000 0.000 0.000;
1.000 1.000 0.077;
1.000 1.000 1.000;
0.500 1.000 1.000;
0.020 1.000 0.117;
1.000 1.000 1.000;
0.020 0.036 1.000;
0.296 0.571 0.113];

w=W(:,refindex)';
S=codonise64(aln.seq);

% Exclude those codons whose w=0 by replacing them with ATG (w=1);
excld=find(w==0);
S(find(ismember(S,excld)))=15;	 % ATG

[n,m]=size(S);
X=w(S);
cai=exp(sum(log(X),2)./m);