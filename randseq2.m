function [s]=randseq2(n,m,seqtype,freq)
%RANDSEQ - Generates a random sequence
%
% Syntax: [s]=randseq2(n,m,seqtype,freq)
%
% Inputs:
%    n         - Number of seq
%    m         - Length of seq
%    seqtype   - 1-nucleotide, 2-codon, 3-protein
%    freq      - Expected frequencies of alphabets
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


if nargin<4
    freq=[.25, .25, .25, .25];
end
if nargin<3
    seqtype=1;
end
if nargin<2
    m=100;
    disp ('100 base generated.')
end

if nargin<1
    n=10;
    disp ('10 sequences generated.')
end

if ~(n>0 && m>0), error('n and m must greater than zero.'); end
%check seqtype and freq combinations

switch (seqtype),
    case (1),
        s = randi([1, 4], n, m);
    case (2),
        %freq_aspergillus_fumigatus=[14.2 23.5 34.11 14.42 11.86 20.98 12.13 13.6 7.37 15.63 5.91 9.52 4.59 28.21 21.08 17.32 14.47 12.9 25.43 11.41 11.29 17.55 12.91 16.36 9.39 15.49 10.83 10.88 6.54 23.43 25.17 15.4 23.09 29.96 36.85 26 14.98 28.87 18.18 22.19 14.32 25.61 10.29 18.7 5.5 25.51 16.84 14.6 18.44 11.44 9.98 18.23 14.9 13.47 7.91 14.71 4.76 4.09 26.13 15.58 13.09];
	%freq=freq_aspergillus_fumigatus./(sum(freq_aspergillus_fumigatus));
	s = decodonise61(randi([1, 61], n, m));
    case (3),
        s = randi([1, 20], n, m);
 end


function [y] = decodonise61(x)
[n,m]=size(x);
y=ones(n,m*3)*5;
[TABLE,CODON] = codontable;
icode=1;
stops=TABLE(icode,:)=='*';
	%CODON=[1 1 1; 1 1 2; 1 1 3; 1 1 4; 1 2 1; 1 2 2; 1 2 3; 1 2 4;
	%       1 3 1; 1 3 2; 1 3 3; 1 3 4; 1 4 1; 1 4 2; 1 4 3; 1 4 4;
	%       2 1 1; 2 1 2; 2 1 3; 2 1 4; 2 2 1; 2 2 2; 2 2 3; 2 2 4;
	%       2 3 1; 2 3 2; 2 3 3; 2 3 4; 2 4 1; 2 4 2; 2 4 3; 2 4 4;
	%       3 1 1; 3 1 2; 3 1 3; 3 1 4; 3 2 1; 3 2 2; 3 2 3; 3 2 4;
	%       3 3 1; 3 3 2; 3 3 3; 3 3 4; 3 4 1; 3 4 2; 3 4 3; 3 4 4;
	%       4 1 1; 4 1 2; 4 1 3; 4 1 4; 4 2 1; 4 2 2; 4 2 3; 4 2 4;
	%       4 3 1; 4 3 2; 4 3 3; 4 3 4; 4 4 1; 4 4 2; 4 4 3; 4 4 4];
	CODON(stops,:)=[];

for i=1:n
for j=1:m
      c=CODON(x(i,j),:);
      y(i,(j-1)*3+1)=c(1,1);
	  y(i,(j-1)*3+2)=c(1,2);
	  y(i,(j-1)*3+3)=c(1,3);
end
end
