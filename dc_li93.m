function [dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln)
%DC_LI93 - Computes Syn. non-syn. substitutions rates by Li93 method
%This is a method for computing the numbers of synonymous and nonsynonymous
%substitutions and the numbers of potentially synonymous and potentially
%nonsynonymous sites based on (Li 1993).
%
% Syntax: [dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln)
%
% Inputs:
%    aln     - Alignment structure
%
% Outputs:
%    dS      - Synonymous substitution rate
%    dN      - Nonsynonymous substitution rate
%    dN_dS   - Ratio of dN and dS
%
% See also: DC_NEI_GOJOBORI86, DC_LI_WU_LUO85, COUNTDEGENERATESITES

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (isstruct(aln)),
	seq=aln.seq;
	geneticcode=aln.geneticcode;
else
	seq=aln;
	geneticcode=1;
end

%if ~(isvalidaln(aln,'CODING')), error ('ERROR: Not protein-coding sequence'); end
if (hasgap(seq))
	disp ('WARNING: Gaps in sequences have been removed.')
	seq=rmcodongaps(seq);
end

% disp('WARNING: This is a simplified implementation of the method.');

[n,m]=size(seq);

L0=zeros(n); L2=zeros(n); L4=zeros(n);
D0=zeros(n); D2=zeros(n); D4=zeros(n);
dS=zeros(n); dN=zeros(n); dN_dS=zeros(n);


for i=1:n-1
for j=i+1:n
%	[seq0,seq2,seq4,L0(i,j),L2(i,j),L4(i,j)]=extractdegeneratesites([aln2.seq(i,:);aln2.seq(j,:)]);
%	fprintf(['L0 %d,L2 %d,L4 %d\n'],L0(i,j),L2(i,j),L4(i,j));
	% [seq0,seq2,seq4,m0,m2,m4]=countdegeneratesites([aln2.seq(i,:);aln2.seq(j,:)]);
	% [seq0,seq2,seq4,m0,m2,m4]=countdegeneratesites(aln2.seq([i j],:));
	[p0,q0,p2,q2,p4,q4,m0,m2,m4]=countdegeneratesites(seq([i j],:));

	% [p0,q0]=countseqpq(seq0); [p2,q2]=countseqpq(seq2); [p4,q4]=countseqpq(seq4);
	% fprintf(['P0 %d %d P2 %d %d P4 %d %d\n'],P0,Q0,P2,Q2,P4,Q4);
	% [x,mx0]=size(seq0); [x,mx2]=size(seq2); [x,mx4]=size(seq4);
	% p0=p0/mx0; q0=q0/mx0; p2=p2/mx2; q2=q2/mx2; p4=p4/mx4; q4=q4/mx4;



	if ((1-2*p0-q0)~=0), a0=1/(1-2*p0-q0); else a0=realmax; end
	if ((1-2*p2-q2)~=0), a2=1/(1-2*p2-q2); else a2=realmax; end
	if ((1-2*p4-q4)~=0), a4=1/(1-2*p4-q4); else a4=realmax; end
	if ((1-2*q0)~=0), b0=1/(1-2*q0); else b0=realmax; end
	if ((1-2*q2)~=0), b2=1/(1-2*q2); else b2=realmax; end
	if ((1-2*q4)~=0), b4=1/(1-2*q4); else b4=realmax; end


	% REF: Nei & Kumar 2000, page 62. (4.9a and 4.9b)
	A0=(1/2)*log(a0)-(1/4)*log(b0); B0=(1/2)*log(b0);
	A2=(1/2)*log(a2)-(1/4)*log(b2); B2=(1/2)*log(b2);
	A4=(1/2)*log(a4)-(1/4)*log(b4); B4=(1/2)*log(b4);

	L0(i,j)=m0/2; L2(i,j)=m2/2; L4(i,j)=m4/2;
	D0(i,j)=A0+B0; D2(i,j)=A2+B2; D4(i,j)=A4+B4;

     	% REF: Nei & Kumar 2000, page 62. (4.10a and 4.10b)
	dS(i,j) = x2real((m2*A2+m4*A4)/(m2+m4)+B4);
	dN(i,j) = x2real((m2*B2+m0*B0)/(m2+m0)+A0);
	dN_dS(i,j)=x2real(dN(i,j)/dS(i,j));

L0(j,i)=L0(i,j); L2(j,i)=L2(i,j); L4(j,i)=L4(i,j);
D0(j,i)=D0(i,j); D2(j,i)=D2(i,j); D4(j,i)=D4(i,j);
dS(j,i)=dS(i,j); dN(j,i)=dN(i,j); dN_dS(j,i)=dN_dS(i,j);
end
end


%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [x2] = x2real(x)
if(~isreal(x))
	x2=nan;
else
	x2=x;
end