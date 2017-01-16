function [v,V] = codonvolatility(aln)
%CODONVOLATILITY - Calculates codon volatility
%REFERENCE: JOSHUA B. PLOTKIN, JONATHAN DUSHOFF & HUNTER B. FRASER, Detecting
%selection using a single genome sequence of M. tuberculosis and P. falciparum
%
% Syntax:  [v,V] = codonvolatility(aln)
%
% Inputs:
%    aln    - Input sequences (Alignment struct)
%
% Outputs:
%    v   - Mean codon volatility for input sequences
%    V   - Matrix containing codon volatility for each codon
%
% See also: GETSYNNONSYNSITES, GETSYNNONSYNDIFF

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if ~(isvalidaln(aln,'CODING'))
	error ('ERROR: Not coding seq')
end
if (hasgap(aln))
	disp ('WARNING: Sequences cannot contain gaps. Gaps will be removed by using REMOVEGAPS.')
	aln=rmcodongaps(aln);
end

global CODON TABLE stops icode

h = waitbar(0.15,'Please wait...');

icode=aln.geneticcode;
[CodonSeq]=codonise64(aln.seq);

[n,m3]=size(aln.seq);
m=m3/3;
V=zeros(n,m);

[TABLE,CODON] = codontable;
stops=find(TABLE(icode,:)=='*');



if (m>64)
	SN=zeros(1,64);
	for (k=1:64), [SN(1,k)] = cal_v(k); end
	V=SN(CodonSeq);
else
	for i=15:30, waitbar(i/100,h); end
	for (i=1:n),
	      seq=CodonSeq(i,:);
		for (j=1:m),
		     cdon=seq(1,j);
		     V(i,j)=V(i,j)+cal_v(cdon);
		end % for (j=1:m)
	end
	for i=30:95, waitbar(i/100,h); end
end

v = sum(V,2)./m;


waitbar(1,h);
close(h);



%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [v] = cal_v(cdon)
%Calculate codon volatility for single codon

global CODON TABLE stops icode

     if ~(ismember(cdon,stops))
	SynDif=0;
	AsynDif=0;
     	for (k=1:64)
		if ~(ismember(k,stops))
			ndiff=sum(CODON(cdon,:)~=CODON(k,:));
			if (ndiff==1)
				SynDif = SynDif + (TABLE(icode,k)==TABLE(icode,cdon));
				AsynDif = AsynDif + (TABLE(icode,k)~=TABLE(icode,cdon));
			end
		end
	end
	v = AsynDif./(AsynDif+SynDif);
     else
	v = 0;
     end