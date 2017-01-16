function [opt,seq1,seq2] = needlemanwunsch(x,y,S,g,fg)
%NEEDLEMANWUNSCH - Calculates optimal gobal alignment of two proteins.
%
% Syntax:  [opt,seq1,seq2] = needlemanwunsch(x,y,S,g,fg)
%
% Inputs:
%    x      - string of amino acids                                   
%    y      - string of amino acids                                   
%    S      - 20 x 20 scoring matrix                                  
%    g      - gap penalty (should be a negative number)               
%    fg     - penalty for the first gap (should be a negative number)
%
% Outputs:
%    opt    - score of the optimal alignment
%    seq1   - x, aligned to y (possibly containing gaps)
%    seq2   - y, aligned to x (possibly containing gaps)
%
% Amino acids corresponding to rows and columns of S must be in the following order:
% "CSTPAGNDEQHRKMILVFYW". Matrix S should be symmetric. The function does NOT test symmetry.
% Case of the characters in p1 and p2 does not influence the final alignment.
% Letters 'B' and 'Z' are treated as 'D' and 'E', respectively. Every amino acid
% scores -1 if aligned with X. If there are more than one alignments with the same
% score, the function returnes only one of them.
%
% Example: p1 = 'ACCALI'
%          p2 = 'ACQVI'
%          S = BLOSUM62
%          g = -4
%          fg = -12
%
%          The function returns:
%
%          d = 5
%          seq1 = ACCALI
%          seq2 = A-CQVI
%
% For information about BLOSUM62 matrix see S. Henikoff, J. Henikoff, "Amino acid
% substitution matrices from protein blocks," Proc. Natl. Acad. Sci. USA, Vol 89.,
% pp. 10915-10919, November 1992, Biochemistry.
%
%
% See also: 

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005


if (nargin<5), g = -4;  end % The gap-open penalty (a nonnegative integer)
if (nargin<4), fg = -12; end % The gap-extend penalty (a nonnegative integer)
if (nargin<3), S = i_BLOSUM62; end
if (nargin<2), error('This function need two sequences.'); end


%S=i_score(x,y,d);

% edit the score matrix
n=length(x); m=length(y);
% D=[zeros(1,m+1); zeros(n,1) S];
D=S;
T=zeros(n,m);

for (i=2:n),
     for (j=2:m),
	[D(i,j),T(i-1,j-1)] = max([0; D(i-1,j-1)+D(i,j);...
	  D(i-1,j)-g; D(i,j-1)-g]);
     end
end



for (i=2:(n+1)),
     for (j=2:(m+1)),
	[D(i,j),T(i-1,j-1)] = max([0; D(i-1,j-1)+D(i,j);...
	  D(i-1,j)-g; D(i,j-1)-g]);
     end
end
D=D(2:end,2:end);

% Trace back the optimal score.
[opt,I]=max(D(1:end));
j=ceil(I/n); opt_y=[];
i=I-n*(j-1); opt_x=[];

nt_x=[]; nt_y=[];
NT=['A','C','G','T','-'];

while (i*j>0&D(i,j)>0),
  opt_y=[j opt_y];
  opt_x=[i opt_x];

  switch (T(i,j))
      case (2)
	  nt_x=[NT(x(i)) nt_x];
	  nt_y=[NT(y(j)) nt_y];
           i=i-1; j=j-1;
       case (3)
	   nt_y=[NT(5) nt_y];
	   nt_x=[NT(x(i)) nt_x];
	   i=i-1;
      case (4)
   	   nt_x=[NT(5) nt_x];
	   nt_y=[NT(y(j)) nt_y];
           j=j-1;
  end
end


% Plot the dot matrix
contour(D)
image(D)
hold on
plot(opt_y, opt_x,'*red')
hold off

disp([nt_x; nt_y])



function [s] = i_score(x,y,d)
% Return score matrix for x & y in {1,2,3,4}
% d = mismatch penalty
%
%Example:
% >>x=ceil(4*rand(1,5))
% >>y=ceil(4*rand(1,7))
% >>[0 y; x' i_score(x,y,1)

%s=(double(x==1)'*double(y==1)+double(x==2)'*double(y==2)+double(x==3)'*double(y==3) ...
%+ double(x==4)'*double(y==4))*(1+d)-d;

s=(double(x==1)'*double(y==1)+double(x==2)'*double(y==2)+double(x==3)'*double(y==3) ...
+ double(x==4)'*double(y==4))*(d);


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
