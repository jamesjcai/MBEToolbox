function [D,VarD]=dn_f84(aln,freq)
%DN_F84 - Distance of Felsenstein 84 model
%
% Syntax: [D,VarD]=dn_f84(aln,freq)
%
% Inputs:
%    aln    - Alignment structure
%    freq   - (optional) 1x4 vector of equilibrium base frequencies
%
% Outputs:
%    D      - Distance matrix
%    VarD   - Variance of distance
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
empirical = 1;

if (nargin==2)
	%[x,y]=size(freq);
	%	if ~(x==1&y==4),
	%		error('not valid freq.')
	%	end
	%	if (sum(freq)-1>eps),
	%		error('not valid freq.')
	%	end
	freq=i_assertfreq(freq);
	% use empirical base frequencies?
	empirical = 0;
end


for (i=1:n),
for (j=i:n),
	if ~(i==j)
	S1=S(i,:); S2=S(j,:);
	if (empirical==1)
		[freq] = estimatefreq(S);
		[d,v,k] = i_SeqPairDistanceF84([S1;S2]);
	else
		[d,v,k] = i_SeqPairDistanceF84([S1;S2],freq);
    end

	D(i,j)=d(1,2); D(j,i)=d(1,2);

    if isnan(v)
	VarD(i,j)=nan; VarD(j,i)=nan;
    else
	VarD(i,j)=v(1,2); VarD(j,i)=v(1,2);
    end

	end
end
end




function [d,v,k] = i_SeqPairDistanceF84(S,freq)
%   This calculates kappa and d from P (proportion of transitions) & Q
%   (proportion of transversions) & pi under F84.
%   When F84 fails, we try to use K80.  When K80 fails, we try
%   to use JC69.  When JC69 fails, we set distance t to maxt.
%   Variance formula under F84 is from Tateno et al. (1994), and briefly
%   checked against simulated data sets.

d=nan; k=nan; v=nan;

[n,m]=size(S);
[P,Q]=countseqpq(S); P=P./m; Q=Q./m;

% Qsmall=min2(1e-10,0.1/m)

if (nargin==2)
	Pi=freq;
else
	N=[sum(sum(S==1)),sum(sum(S==2)),sum(sum(S==3)),sum(sum(S==4))];
	Pi=N./sum(N);		% Pi(A,C,G,T);
end


tc=Pi(4)*Pi(2);
ag=Pi(1)*Pi(3);
R=Pi(1)+Pi(3);
Y=Pi(2)+Pi(4);
A=tc/Y+ag/R; B=tc+ag; C=Y*R;

a=(2*B+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*C)) - P) / (2*A);
b=1-Q/(2*C);
%a=a(1,2);
%b=b(1,2);

if (any(a(:)<=0) || any(b(:)<=0))
  failF84=1;
  return;
end

a=-.5*i_safelog(a); b=-.5*i_safelog(b);
if(b<=0)
  failF84=1;
	  % try to use K80
	  W1 = 1-2*P-Q; W2 = 1-2*Q;
	  d=(-1/2)*i_safelog(W1)-(1/4)*i_safelog(W2);
 	  if (nargout==2)
		W1=1./W1; W2=1./W2; W3=(W1+W2)./2;
		v=((W1.^2).*P + (W3.^2).*Q - (W1.*P+W3.*Q).^2)./m;
	  end
	  %end of 'try to use K80'
  return;
end

k_F84=a/b-1;

if (k_F84>999)
	k_F84=999;
end
k_F84=max(k_F84,-0.5);

k=k_F84;
d = 4*b*(tc*(1+ k_F84/Y)+ag*(1+ k_F84/R)+C);

if (nargout>1)
     a = A.*C./(A.*C-C.*P/2-(A-B).*Q./2);
     b = A.*(A-B)./(A.*C-C.*P./2-(A-B).*Q./2)- (A-B-C)./(C-Q./2);
     v = (a.*a.*P+b.*b.*Q-(a.*P+b.*Q)^2)./n;
end