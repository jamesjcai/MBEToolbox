function [SynDif,AsynDif] = getsynnonsyndiff(icode);
%GETSYNNONSYNDIFF - Return matrices of Syn- Nonsyn- differences between codons
%
% Syntax: [SynDif,AsynDif] = getsynnonsyndiff(icode)
%         [SynDif,AsynDif] = getsynnonsyndiff;
%
% Inputs:
%    icode      - index of genetic code
%
% Outputs:
%    SynDif      - 64x64 matrix of syn- differences between codons
%    AsynDif     - 64x64 matrix of nonsyn- differences between codons
%
%
% See also: GETSYNNONSYNSITES, DC_NG86

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if(nargin<1), icode=1; end
if (icode==1)
	load synnonsyndiff SynDif AsynDif;
	return;
end

SynDif=zeros(64,64);
AsynDif=zeros(64,64);
[TABLE,CODON] = codontable;


% load SNMatrix SN ns na;
% index=1:64;
% stops=index(logical(TABLE(icode,:)=='*'));

stops=find(TABLE(icode,:)=='*');

for (i=1:64)
for (j=i:64)
		[SynDif(i,j),AsynDif(i,j)] = i_countSynNonsynDiff(i,j,stops,...
		                             CODON,icode,TABLE);
		SynDif(j,i)=SynDif(i,j); AsynDif(j,i)=AsynDif(i,j);
end
end


%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

function [SynDif,AsynDif] = i_countSynNonsynDiff(codon1,codon2,stops,CODON,...
                            icode,TABLE)

index=1:64;
SynDif=0; AsynDif=0;

if (ismember(codon1,stops)|ismember(codon2,stops))
	SynDif=-1; AsynDif=-1;
	return;
end


ndiff=sum(CODON(codon1,:)~=CODON(codon2,:));
switch (ndiff)
    case (0)
	 SynDif=0; AsynDif=0;
    case (1)
	SynDif = TABLE(icode,codon1)==TABLE(icode,codon2);
	AsynDif = TABLE(icode,codon1)~=TABLE(icode,codon2);
    case (2)
	A1=CODON~=repmat(CODON(codon1,:),64,1);
	A2=CODON~=repmat(CODON(codon2,:),64,1);
	A12 = intersect(find(sum(A1')==1),find(sum(A2')==1));
	A13 = intersect(find(sum(A1')==1),find(sum(A2')==1));

	A12 = setdiff(A12,stops);
	[n,m]=size(A12);

	% if(m~=2) then m should = 1, disp('stop codon in pathway!'); end
	if (m==0)
		SynDif=0; AsynDif=2;
		return;
	end


	SynDif = sum(TABLE(icode,repmat(codon1,1,m))==TABLE(icode,A12)) +...
		 sum(TABLE(icode,repmat(codon2,1,m))==TABLE(icode,A12));
	AsynDif = sum(TABLE(icode,repmat(codon1,1,m))~=TABLE(icode,A12)) +...
		  sum(TABLE(icode,repmat(codon2,1,m))~=TABLE(icode,A12));
	SynDif=SynDif/m;
	AsynDif=AsynDif/m;
    case (3)

	A1=CODON~=repmat(CODON(codon1,:),64,1);
	A2=CODON~=repmat(CODON(codon2,:),64,1);
	% X1=index((sum(A1')==1)'); X2=index((sum(A2')==1)');
	% Y1=index((sum(A1')==2)'); Y2=index((sum(A2')==2)');

	X1=find(sum(A1')==1); X2=find(sum(A2')==1);
	Y1=find(sum(A1')==2); Y2=find(sum(A2')==2);
	Z1 = intersect(X1,Y2);	Z2 = intersect(X2,Y1);
	% Z1 = setdiff(Z1,stops);	Z2 = setdiff(Z2,stops);
	Z1=Z1([1 1 2 2 3 3]);		Z2=Z2([1 2 1 3 2 3]);
	codon1_6 = codon1([1 1 1 1 1 1]); 	codon2_6 = codon2([1 1 1 1 1 1]);
	% SixPath = [codon1;Z1;Z2;codon2];
	SixPath = cat(1,codon1_6,Z1,Z2,codon2_6);
	picker=sum(~ismember(SixPath,stops))==4;	% pathway contains four codon, ie, without stop codons

	SixPath=SixPath(:,picker);
	[n,m]=size(SixPath);
	if (m==0)
		% codon1
		% codon2
		% CODON(codon1,:)
		% CODON(codon2,:)
		% sum(CODON(codon1,:)~=CODON(codon2,:))
		% fprintf(['%d %d\n'], codon1,codon2);
		SynDif=1; AsynDif=2;
		return;
	end
	SynDif = sum(TABLE(icode,SixPath(1,:))==TABLE(icode,SixPath(2,:))) +...
		 sum(TABLE(icode,SixPath(2,:))==TABLE(icode,SixPath(3,:))) +...
		 sum(TABLE(icode,SixPath(3,:))==TABLE(icode,SixPath(4,:)));
	AsynDif = sum(TABLE(icode,SixPath(1,:))~=TABLE(icode,SixPath(2,:))) +...
		  sum(TABLE(icode,SixPath(2,:))~=TABLE(icode,SixPath(3,:))) +...
		  sum(TABLE(icode,SixPath(3,:))~=TABLE(icode,SixPath(4,:)));
	SynDif=SynDif/m;
	AsynDif=AsynDif/m;
end