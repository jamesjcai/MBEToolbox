function [DS,DN]=normaliseKaKs(Aln,nboot)

if nargin < 2 | isempty(nboot) | isnan(nboot), nboot = 1000; end

if ~(isValidAln(Aln,'CODING')), error ('ERROR: Not protein-coding sequence'); end
if (hasGap(Aln))
	% disp ('WARNING: Gaps in sequences have been removed.')
	Aln=removeGappedCodon(Aln);
end


% Aln2 = copyAlnHeader(Aln);
[n,m3]=size(Aln.seq);
m=m3/3;
[CodonSeq]=codoniseSeq(Aln.seq);

[ks1,ka1]=i_dc_nei_gojobori86(CodonSeq);

if ~(n==2)
	error('Must be two sequences.')
end


DS=ones(1,nboot)*999;
DN=ones(1,nboot)*999;
counter=0;
sumks=0;
sumka=0;

 h = waitbar(0,'Please wait...');
%for (i=1:nboot),
while (counter<nboot),
	index = randperm(m);
        C1= CodonSeq(1,index);
	C2=CodonSeq(2,:);
	CS = [C1;C2];
	
	% CS=realign(CS)
	[ks,ka]=i_dc_nei_gojobori86(CS);
	if (ks>0&ks<1)
		counter=counter+1;
    		DS(counter)=ks;
		DN(counter)=ka;
		sumks=sumks+ks;
		sumka=sumka+ka;
	end

 waitbar(counter/nboot,h);
end

 waitbar(1,h);
 close(h);


DS(find(DS==999))=[];
DN(find(DN==999))=[];

if (counter<3)
    ds=999;
    dn=999;
    Zs=999;
    Zn=999;     
else 
    ds=sumks/counter;
    dn=sumka/counter;
    
if std(DS)==0
    Zs=999;
else    
    Zs = ks1-mean(DS)/std(DS);
end
if std(DN)==0
    Zn=999;
else    
    Zn = ka1-mean(DN)/std(DN);
end
    
    
end
% n=counter;



fprintf(['%4.6f	%4.6f	%4.6f	%4.6f\n'],ds,dn,Zs,Zn);




function [ds,dn,dS,dN]=i_dc_nei_gojobori86(CodonSeq)

[n,ignore] = size(CodonSeq);
global ns na

[S,N] = getSynNonsynSites(1);
[ns,na] = getSynNonsynDiff(1);

dS=zeros(n); dN=zeros(n); dN_dS=zeros(n);
St=zeros(n); Nt=zeros(n); nst=zeros(n); nat=zeros(n);
VdS=zeros(n);VdN=zeros(n);

i=1;
j=2;

	if (i~=j)
		St(i,j) = sum(sum(S(CodonSeq([i j],:))))/2;
		Nt(i,j) = sum(sum(N(CodonSeq([i j],:))))/2;
		[nst(i,j), nat(i,j)] = calculate_diff(CodonSeq([i j],:));
		pS=nst(i,j)/St(i,j); 
		pN=nat(i,j)/Nt(i,j);

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

	ds=dS(i,j); dn=dN(i,j);



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