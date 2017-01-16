function [D,VarD,K]=dn_khy(aln,freq)
%DN_HKY - Hasegawa, Kishino and Yano 85 (HKY) Distance
%
% Syntax: [D,VarD]=dn_khy(aln,freq)
%
% Inputs:
%    aln    - Alignment structure
%    freq   - (optional) 1x4 vector of equilibrium base frequencies
%
% Outputs:
%    D      - Distance matrix
%    VarD   - Variance of distance
%    K      - kappa
%
%REF: Yan Z and Nielsen R (2000) Estimating synonymous and nonsynonymous
%     subsititution rates under realistic evolutionary models. p. 34
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


if (isstruct(aln)),
	S=aln.seq;
else
	S=aln;
end

[n,m] = size(S);
if (n<2)
	error('At least two sequences.')
end

D=zeros(n,n);
VarD=zeros(n,n);
K=zeros(n,n);
empirical = 0;

if (nargin==2)
	[x,y]=size(freq);
	if ~(x==1&&y==4),
		error('not valid freq.')
	end
	if (sum(freq)-1>eps),
		error('not valid freq.')
	end
	% use empirical base frequencies?
	empirical = 1;
end


for (i=1:n),
for (j=i:n),
	if ~(i==j)
	S1=S(i,:); S2=S(j,:);
	if (empirical==1)
		[d,v,k] = i_SeqPairDistanceHKY85([S1;S2],freq);
	else
		[freq] = estimatefreq(S);
		[d,v,k] = i_SeqPairDistanceHKY85([S1;S2],freq);
	end

    K(i,j)=k; K(j,i)=k;
	D(i,j)=d; D(j,i)=d;
	VarD(i,j)=v; VarD(j,i)=v;
	end
end
end




function [d,v,k] = i_SeqPairDistanceHKY85(S,freq)
%   This calculates kappa and d from P (proportion of transitions) & Q
%   (proportion of transversions) & pi under F84.
%   When F84 fails, we try to use K80.  When K80 fails, we try
%   to use JC69.  When JC69 fails, we set distance t to maxt.
%   Variance formula under F84 is from Tateno et al. (1994), and briefly
%   checked against simulated data sets.

d=nan; k=nan; v=nan;

[n,m]=size(S);
[P1,P2,Q]=countseqppq(S);
P1=P1./m; P2=P2./m; Q=Q./m;

% Qsmall=min2(1e-10,0.1/m)

if (nargin==2)
	Pi=freq;
else
	N=[sum(sum(S==1)),sum(sum(S==2)),sum(sum(S==3)),sum(sum(S==4))];
	Pi=N./sum(N);		% Pi(A,C,G,T);
end


tc=Pi(4)*Pi(2);
ag=Pi(1)*Pi(3);

%   Y=p[0]+p[1];    R=p[2]+p[3];  tc=p[0]*p[1]; ag=p[2]*p[3];
R=Pi(1)+Pi(3);
Y=Pi(2)+Pi(4);
A=tc/Y+ag/R; B=tc+ag; C=Y*R;

% HKY85 model, following Tamura (1993, MBE, ..)
%   alpha=0 if no gamma
%   return -1 if in error.
%   Check DistanceF84() if variances are wanted.
% 0=t 1=c 2=a 3=g
%   return 4*p[0]*p[1]*a1 + 4*p[2]*p[3]*a2 + 4*Y*R*b;
%   return 4*p[T]*p[C]*a1 + 4*p[A]*p[G]*a2 + 4*Y*R*b;



      a1=1-Y*P1/(2*tc)-Q/(2*Y);  a2=1-R*P2/(2*ag)-Q/(2*R);   b=1-Q/(2*Y*R);
      if (a1<=0 || a2<=0 || b<=0)
	   d=nan;       return;
      end
      a1=-log(a1); a2=-log(a2); b=-log(b);
      a1=.5/Y*(a1-R*b);  a2=.5/R*(a2-Y*b);  b=.5*b;
      d= 4*Pi(4)*Pi(2)*a1 + 4*Pi(1)*Pi(3)*a2 + 4*Y*R*b;

if (nargout>1)
     P=P1+P2;
     a = A*C/(A*C-C*P/2-(A-B)*Q/2);
     x = A*(A-B)/(A*C-C*P/2-(A-B)*Q/2) - (A-B-C)/(C-Q/2);
     v = (a*a*P+x*x*Q-(a*P+x*Q)^2)/n;
end

if (nargout>2)
      k=999;
      if (b>0)
          k = min((a1+a2)/2/b, 999);
      end
end




