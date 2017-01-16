function [S,N] = getsynnonsynsites_weighted(codonusenum,icode);
%GETSYNNONSYNSITES - Return matrix of Syn- Nonsyn- sites of codons
%
% Syntax: [S,N] = getsynnonsynsites(icode)
%         [S,N] = getsynnonsynsites;
%
% Inputs:
%    icode      - index of genetic code
%    codonusenum - 1x64 vector of codon usage number
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



if(nargin<2), icode=1; end
if(nargin<1), codonusenum=ones(1,64); end



% global CODON TABLE stops
S=zeros(1,64); N=zeros(1,64);
[ct,CODON] = codontable(icode);
stops=find(ct=='*');


for (k=1:64),
	[S(1,k),N(1,k)] =...
        i_countSynNonsynSite(k,ct,CODON,stops,codonusenum);
end



%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%
function [SynSite,AsynSite] = ...
    i_countSynNonsynSite(codon1,ct,CODON,stops,codonusenum)
    % global CODON TABLE stops

SynSite=0; AsynSite=3;
if (ismember(codon1,stops)) return; end

% repmat(CODON(codon1,:),64,1) constructs a repeated copies of input codon
% A3 marks different base of 64 codons from the input codon
A3=CODON~=repmat(CODON(codon1,:),64,1);

neighbor1=find(sum(A3,2)==1&A3(:,1)==1);
neighbor2=find(sum(A3,2)==1&A3(:,2)==1);
neighbor3=find(sum(A3,2)==1&A3(:,3)==1);

nonstop_neighbor1=setdiff(neighbor1,stops);
nonstop_neighbor2=setdiff(neighbor2,stops);
nonstop_neighbor3=setdiff(neighbor3,stops);


rc=codonusenum(nonstop_neighbor1)./sum(codonusenum(nonstop_neighbor1));
if any(ct(codon1)==ct(nonstop_neighbor1))
    idx=find(ct(codon1)==ct(nonstop_neighbor1));
    SynSite1=sum(rc(idx));
else
    SynSite1=0;
end


rc=codonusenum(nonstop_neighbor2)./sum(codonusenum(nonstop_neighbor2));
if any(ct(codon1)==ct(nonstop_neighbor2))
    idx=find(ct(codon1)==ct(nonstop_neighbor2));
    SynSite2=sum(rc(idx));
else
    SynSite2=0;
end

rc=codonusenum(nonstop_neighbor3)./sum(codonusenum(nonstop_neighbor3));
if any(ct(codon1)==ct(nonstop_neighbor3))
    idx=find(ct(codon1)==ct(nonstop_neighbor3));
    SynSite3=sum(rc(idx));
else
    SynSite3=0;
end

%SynSite3=sum(ct(codon1)==ct(nonstop_neighbor3))/length(nonstop_neighbor3);


SynSite=SynSite1+SynSite2+SynSite3;
AsynSite=3-SynSite;
