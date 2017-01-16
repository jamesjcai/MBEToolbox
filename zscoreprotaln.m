function [z,p1,p2,score_ori,vscore1,vscore2] = zscoreprotaln(filename,nboot)
%ZSCOREPROTALN - Z score of protein alignment
%The Z-score gives a measure of significance against a background of randomly
%generated sequences with the same composition and length as the original
%sequences. Hence it is designed to overcome the bias due to the composition
%of the alignment. A Z-score is usually calculated by comparing an actual
%alignment score with the scores obtained on a set of random sequences by a
%Monte-Carlo process. The Z-score is defined as: $$Z(A,B)=(S(A,B)-mean)/standard
%\mbox{ } deviation$$ where $S(A,B)$ is the Smith-Waterman (S-W) score between
%wo sequence $A$ and $B$, the mean and standard deviation are taken with realigning
%the permuted sequences.
%
%@article{
%3162770,
%   Author = {Pearson, W. R. and Lipman, D. J.},
%   Title = {Improved tools for biological sequence comparison},
%   Journal = {Proc Natl Acad Sci U S A},
%   Volume = {85},
%   Number = {8},
%   Pages = {2444-8},
%   Accession = {3162770},
%   Year = {1988} }
%
% Syntax: [z,p1,p2,score_ori,vscore1,vscore2] = zscore_protalign(filename,nboot)
%
% Inputs:
%    filename  - file
%    nboot     - number of bootstrapping
%
% Outputs:
%    z           - Z score
%    p1          - protein sequence 1 string
%    p2          - protein sequence 2 string
%    score_ori   - Orignial raw score
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



% g = gap penalty (should be a negative number)
% fg = penalty for the first gap (should be a negative number)

g = -4;  % The gap-open penalty (a nonnegative integer in map.c)
fg = -12; % The gap-extend penalty (a nonnegative integer in map.c)
S = i_BLOSUM62;
[p1,p2] = i_readFASTA(filename);

m1=length(p1);
m2=length(p2);

score_ori = protglobal (p1, p2, S, g, fg);
vscore1=zeros(1,nboot);
vscore2=zeros(1,nboot);

for (k=1:nboot),
	p1_rep = p1(:,randperm(m1));
	vscore1(1,k) = protglobal (p1_rep, p2, S, g, fg);
end

for (k=1:nboot),
	p2_rep = p2(:,randperm(m2));
	vscore2(1,k) = protglobal (p1, p2_rep, S, g, fg);
end

z1=(score_ori-mean(vscore1))/std(vscore1);
z2=(score_ori-mean(vscore2))/std(vscore2);
z=round(min(z1,z2));



function [BLOSUM62] = i_BLOSUM62()
BLOSUM62=[9 -1 -1 -3 0 -3 -3 -3 -4 -3 -3 -3 -3 -1 -1 -1 -1 -2 -2 -2;
-1 4 1 -1 1 0 1 0 0 0 -1 -1 0 -1 -2 -2 -2 -2 -2 -3;
-1 1 5 -1 0 -2 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 0 -2 -2 -2;
-3 -1 -1 7 -1 -2 -2 -1 -1 -1 -2 -2 -1 -2 -3 -3 -2 -4 -3 -4;
0 1 0 -1 4 0 -2 -2 -1 -1 -2 -1 -1 -1 -1 -1 0 -2 -2 -3;
-3 0 -2 -2 0 6 0 -1 -2 -2 -2 -2 -2 -3 -4 -4 -3 -3 -3 -2;
-3 1 0 -2 -2 0 6 1 0 0 1 0 0 -2 -3 -3 -3 -3 -2 -4;
-3 0 -1 -1 -2 -1 1 6 2 0 -1 -2 -1 -3 -3 -4 -3 -3 -3 -4;
-4 0 -1 -1 -1 -2 0 2 5 2 0 0 1 -2 -3 -3 -2 -3 -2 -3;
-3 0 -1 -1 -1 -2 0 0 2 5 0 1 1 0 -3 -2 -2 -3 -1 -2;
-3 -1 -2 -2 -2 -2 1 -1 0 0 8 0 -1 -2 -3 -3 -3 -1 2 -2;
-3 -1 -1 -2 -1 -2 0 -2 0 1 0 5 2 -1 -3 -2 -3 -3 -2 -3;
-3 0 -1 -1 -1 -2 0 -1 1 1 -1 2 5 -1 -3 -2 -2 -3 -2 -3;
-1 -1 -1 -2 -1 -3 -2 -3 -2 0 -2 -1 -1 5 1 2 1 0 -1 -1;
-1 -2 -1 -3 -1 -4 -3 -3 -3 -3 -3 -3 -3 1 4 2 3 0 -1 -3;
-1 -2 -1 -3 -1 -4 -3 -4 -3 -2 -3 -2 -2 2 2 4 1 0 -1 -2;
-1 -2  0 -2 0 -3 -3 -3 -2 -2 -3 -3 -2 1 3 1 4 -1 -1 -3;
-2 -2 -2 -4 -2 -3 -3 -3 -3 -3 -1 -3 -3 0 0 0 -1 6 3 1;
-2 -2 -2 -3 -2 -3 -2 -3 -2 -1 2 -2 -2 -1 -1 -1 -1 3 7 2;
-2 -3 -2 -4 -3 -2 -4 -4 -3 -2 -2 -3 -3 -1 -3 -2 -3 1 2 11];


function [p1,p2] = i_readFASTA(filename)

disp(['Reading ',filename]);

txtc = textread(filename,'%s','delimiter','\n','whitespace','','bufsize', 40950);
mt = find(cellfun('isempty',txtc));
% eliminate empty lines
txtc(mt) = [];

txt=char(txtc);
mt=find(txt(:,1)=='>');
%names=txt(mt,[2:end]);

[n,m]=size(txt);

% add one line at the end of file
mt2=cat(1,mt,n+1);
S=[];
for (k=1:length(mt)),
	sblock = txt([mt2(k)+1:mt2(k+1)-1],:);
	seq = i_mat2vector(sblock);
	S=strvcat(S,upper(seq));
end
p1=removeblanks(S(1,:));
p2=removeblanks(S(2,:));



function [s] = i_mat2vector(S)
s=[];
[n,m]=size(S);
for (k=1:n),
      s=cat(2,s,removeblanks(S(k,:)));
end


function [out] = removeblanks(in)
%REMOVEBLANKS - Removes space within sequence,as well as, both leading and
%trailing blanks
[n,m] = size(in);
codes = zeros(1,m);
for i = 1:m
	for j = 1:n
		if (in(j,i) == ' ')
			codes(i) = 1;
			break;
		end
	end
end
for i = 1:m
	j = m - i + 1;
   if codes(j) == 1
	in(:,j) = [];
   end
end

[r,c] = find( (in~=0) & ~isspace(in) );
if isempty(c),
	out = in([]);
else
	out = in(:,min(c):max(c));
end