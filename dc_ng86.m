function [dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln)
%DC_NG86 - Compute syn. non-syn. substitutions rates by Nei-Gojobori method
%This is a method for computing the numbers of synonymous and nonsynonymous
%substitutions and the numbers of potentially synonymous and potentially
%nonsynonymous sites based on (Nei and Gojobori 1986).
%
% Syntax: [dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln)
%
% Inputs:
%    aln     - Alignment structure
%
% Outputs:
%    dS      - Synonymous substitution rate
%    dN      - Nonsynonymous substitution rate
%    dN_dS   - Ratio of dN and dS
%    VdS     - Variance of dS
%    VdN     - Variance of dN
%    St      - Number of Syn sites
%    Nt      - Number of non-syn sites
%    nst     - Number of Syn differences only
%    nat     - Number of non-Syn differences only
%
% See also: GETSYNNONSYNSITES, GETSYNNONSYNDIFF, DC_LI_WU_LUO85, DC_PAMILO_BIANCHI_LI93

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

[n,m] = size(seq);
[S,N] = getsynnonsynsites(geneticcode);
[ns,na] = getsynnonsyndiff(geneticcode);


dS=zeros(n); dN=zeros(n); dN_dS=zeros(n);
St=zeros(n); Nt=zeros(n); nst=zeros(n); nat=zeros(n);
VdS=zeros(n);VdN=zeros(n);

[CodonSeq]=codonise64(seq);

for i=1:n-1
for j=i+1:n
        xS=S(CodonSeq([i j],:));
        xN=N(CodonSeq([i j],:));
		St(i,j) = sum(sum(xS))/2;
		Nt(i,j) = sum(sum(xN))/2;
		[nst(i,j), nat(i,j)] = calculate_diff(CodonSeq([i j],:),ns,na);
		pS=nst(i,j)/St(i,j); pN=nat(i,j)/Nt(i,j);

		% The approximate large-sample variances of dS and dN can be
		% computed by Equation (3.9) if we replace p in the equation by
		% pS and pN and n by S and N (Nei 1987; Nei and Kumar 2000, p54)
		VdS(i,j)=(pS*(1-pS))/((1-(4/3)*pS)^2*St(i,j));
		VdN(i,j)=(pN*(1-pN))/((1-(4/3)*pN)^2*Nt(i,j));

		%VdS(i,j)=(9*pS*(1-pS))/((3-4*pS)^2*St(i,j));
		%VdN(i,j)=(9*pN*(1-pN))/((3-4*pN)^2*Nt(i,j));

	        dS(i,j) = 1-(4/3)*pS;	% <=== dS
		dN(i,j) = 1-(4/3)*pN;	% <=== dN

			 if(dS(i,j)==1)
				dS(i,j)=0;
			 else
				if (dS(i,j)<0)
					dS(i,j)=-1;
				else
					dS(i,j)=-0.75*log(dS(i,j));
				end
			 end

			 if(dN(i,j)==1)
				dN(i,j)=0;
			 else
				if (dN(i,j)<0)
					dN(i,j)=-1;
				else
					dN(i,j)=-0.75*log(dN(i,j));
				end
			 end

			 if (dS(i,j)<=0)
				dN_dS(i,j) = nan;
			 else
				dN_dS(i,j) = dN(i,j)/dS(i,j);
			 end

	St(j,i)=St(i,j); Nt(j,i)=Nt(i,j);
	nst(j,i)=nst(i,j); nat(j,i)=nat(i,j);
	VdS(j,i)=VdS(i,j); VdN(j,i)=VdN(i,j);
	dS(j,i)=dS(i,j); dN(j,i)=dN(i,j); dN_dS(j,i)=dN_dS(i,j);
end
end

%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [nst, nat] = calculate_diff(CodonSeq,ns,na)
	[n,m]=size(CodonSeq);
	if (n~=2) error ('must be two sequences!'); end
	nst=0;
	nat=0;

	for (k=1:m)
	    xns=ns(CodonSeq(1,k), CodonSeq(2,k));
	    if (xns>0)
		nst = nst+xns;
	    end
	    xna = na(CodonSeq(1,k), CodonSeq(2,k));
	    if (xna>0)
		nat = nat+xna;
	    end
    end

