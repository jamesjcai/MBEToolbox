function [chrid] = ensembl_chrid_v45(gid)
%CHRNUC_V52 - returns nucleotides of reference geneome
%USAGE: chrnuc_v52(chrid,startn,endn)

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007

% 'ENST00000380612'
%if nargin<2, species='human'; end
species='human';
switch lower(species)
    case 'human'
        sptxt='Homo_sapiens';
    case 'mouse'        
        sptxt='Mus_musculus';
    otherwise
        error('Species is not supported.')
end



chrid='';
try
[chrid]=i_getchrid(gid,sptxt);
catch
end



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

%http://www.ensembl.org/Homo_sapiens/Transcript/Export?db=core;g=ENSG00000184809;output=fasta;r=21:39890949-39906238;strand=feature;t=ENST00000380612;time=1235329353.4944;st=coding;_format=Text

% http://jun2007.archive.ensembl.org/Homo_sapiens/exportview?seq_region_name=&type1=transcript&anchor1=ENST00000380612&type2=bp&anchor2=&downstream=&upstream=&format=fasta&action=export&_format=Text&options=coding&options=&output=txt&submit=Continue+%3E%3E

function [chrid] = i_getchrid(gid,sptxt)

%	url=sprintf('http://jun2007.archive.ensembl.org/%s/geneview?gene=%s',...
url=sprintf('http://jun2007.archive.ensembl.org/%s/exportview?_format=Text&seq_region_name=&type1=gene&anchor1=%s&type2=bp&anchor2=&downstream=&upstream=&format=gff&action=export&options=gene&options=&output=txt&submit=Continue+%3E%3E',...
	         sptxt,gid);
       
	[s0,status0]=urlread(url);
	if ~(status0>0),
		disp('Unable to download data.');
		return;
	end
	[s] = strread(s0,'%s','delimiter','\n');
	%s(1)=[];
	%s=char(s);
        
    x=s{1};
    [a,b]=regexp(x,'chromosome:NCBI36:[\dXY]+:');
    chrid=x(1+18:b-1);
    
    %chrid=s{1};
   % strs='';
   % for k=1:size(s,1)
   %     strs=[strs,s(k,:)];
   % end
   % strs=deblank(strs);
    

