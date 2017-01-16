function [S,N] = getsynnonsynsites(icode);
%GETSYNNONSYNSITES - Return matrix of Syn- Nonsyn- sites of codons
%
% Syntax: [S,N] = getsynnonsynsites(icode)
%         [S,N] = getsynnonsynsites;
%
% Inputs:
%    icode      - index of genetic code
%
% Outputs:
%    S,N      - 2x64 matrix of syn- and nonsyn- sites of codons
%
%
% See also: GETSYNNONSYNDIFF, DC_NG86

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if(nargin<1) icode=1; end

if (icode==1)
	load synnonsynsites S N
	return;
end

% global CODON TABLE stops
S=zeros(1,64); N=zeros(1,64);
[TABLE,CODON] = codontable;

stops=find(TABLE(icode,:)=='*');
for (k=1:64),
	[S(1,k),N(1,k)] = i_countSynNonsynSite(k,icode,CODON,TABLE,stops);
end



%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%
function [SynSite,AsynSite] = i_countSynNonsynSite(codon1,icode,CODON,TABLE,stops)
% global CODON TABLE stops

index=1:64;
SynSite=0; AsynSite=3;
if (ismember(codon1,stops)) return; end

A3=CODON~=repmat(CODON(codon1,:),64,1);

% neighbor=index((sum(A3')==1))

neighbor1=index((sum(A3')==1)'&A3(index,1)==1);
neighbor2=index((sum(A3')==1)'&A3(index,2)==1);
neighbor3=index((sum(A3')==1)'&A3(index,3)==1);

nonstop_neighbor1=setdiff(neighbor1,stops);
nonstop_neighbor2=setdiff(neighbor2,stops);
nonstop_neighbor3=setdiff(neighbor3,stops);

[n1,m1]=size(nonstop_neighbor1);
refer1=repmat(codon1,1,m1);
i1=1:m1;
SynSite1=sum(TABLE(icode,refer1(1,i1))==TABLE(icode,nonstop_neighbor1(1,i1)))/m1;


[n2,m2]=size(nonstop_neighbor2);
refer2=repmat(codon1,1,m2);
i2=1:m2;
SynSite2=sum(TABLE(icode,refer2(1,i2))==TABLE(icode,nonstop_neighbor2(1,i2)))/m2;


[n3,m3]=size(nonstop_neighbor3);
refer3=repmat(codon1,1,m3);
i3=1:m3;
SynSite3=sum(TABLE(icode,refer3(1,i3))==TABLE(icode,nonstop_neighbor3(1,i3)))/m3;

SynSite=SynSite1+SynSite2+SynSite3;
AsynSite=3-SynSite;