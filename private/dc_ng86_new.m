function [dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86_new(Aln,ratio)

%Zhang J, Rosenberg HF, Nei M (1998)  Positive Darwinian selection after
%gene duplication in primate ribonuclease genes. Proc. Natl. Acad. Sci.
%USA 95:3708-3713.

%NG-NEW is designed for estimating synonymous and nonsynonymous distances
%between protein coding DNA sequences.  The method is modified from the original 
%Nei and Gojobori (1986) method to take into account the transition bias. 


%DC_NEI_GOJOBORI86 - Compute syn. non-syn. substitutions rates by Nei-Gojobori method
%This is a method for computing the numbers of synonymous and nonsynonymous
%substitutions and the numbers of potentially synonymous and potentially
%nonsynonymous sites based on (Nei and Gojobori 1986).
%
% Syntax: [dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_nei_gojobori86(Aln)
%
% Inputs:
%    Aln     - Alignment structure
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

% Molecular Biology & Evolution Toolbox, (C) 2004
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 8/5/2004


if(ratio<0.000001),
	ratio=1;
end 




if ~(isValidAln(Aln,'CODING')), error ('ERROR: Not protein-coding sequence'); end
if (hasGap(Aln))
	disp ('WARNING: Gaps in sequences have been removed.')
	Aln=removeGappedCodon(Aln);
end

Aln2 = Aln;
[n,m] = size(Aln2.seq);


global ns na

[S,N] = getSynNonsynSites(Aln.geneticcode);
[ns,na] = getSynNonsynDiff(Aln.geneticcode);


dS=zeros(n); dN=zeros(n); dN_dS=zeros(n);
St=zeros(n); Nt=zeros(n); nst=zeros(n); nat=zeros(n);
VdS=zeros(n);VdN=zeros(n);

[CodonSeq]=codoniseSeq(Aln2.seq);

for i=1:n,
for j=i:n,
	if (i~=j)
		St(i,j) = sum(sum(S(CodonSeq([i j],:))))/2;
		Nt(i,j) = sum(sum(N(CodonSeq([i j],:))))/2;
		[nst(i,j), nat(i,j)] = calculate_diff(CodonSeq([i j],:));
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

function [nst, nat] = calculate_diff(CodonSeq)
global ns na
[n,m]=size(CodonSeq);
if (n~=2) error ('must be two sequences!'); end
nst=0;
nat=0;

for (k=1:m)
	nst = nst+ns(CodonSeq(1,k), CodonSeq(2,k));
	nat = nat+na(CodonSeq(1,k), CodonSeq(2,k));
end