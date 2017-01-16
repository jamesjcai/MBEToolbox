function [raws,s] = ensembl_cds(transid,species)
%CHRNUC_V52 - returns nucleotides of reference geneome
%USAGE: chrnuc_v52(chrid,startn,endn)

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007

% 'ENST00000380612'
if nargin<2, species='human'; end

switch lower(species)
    case 'human'
        sptxt='Homo_sapiens';
    case 'chimp'
        sptxt='Pan_troglodytes';
    case 'macaque'        
        sptxt='Macaca_mulatta';        
    case 'mouse'        
        sptxt='Mus_musculus';
    case 'rat'
        sptxt='Rattus_norvegicus';
    case 'dog'
        sptxt='Canis_familiaris';
    case 'zebrafish'            
        sptxt='Danio_rerio';
    otherwise
        error('Species is not supported.')
end
sptxt
[s,raws]=i_getnuc(transid,sptxt);



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

%http://www.ensembl.org/Homo_sapiens/Transcript/Export?db=core;g=ENSG00000184809;output=fasta;r=21:39890949-39906238;strand=feature;t=ENST00000380612;time=1235329353.4944;st=coding;_format=Text

function [strs,s] = i_getnuc(transid,sptxt)

%	base=sprintf('http://www.ensembl.org/%s/Transcript/Export?',sptxt);	
%	base=sprintf('http://aug2006.archive.ensembl.org/%s/exportview?',sptxt);	    
    
    %url=sprintf('http://www.ensembl.org/Homo_sapiens/Transcript/Export?db=core;output=fasta;strand=feature;t=%s;st=coding;_format=Text',...                
    %url=sprintf('http://www.ensembl.org/Homo_sapiens/Transcript/Export?db=core;g=ENSG00000198938;output=fasta;strand=feature;t=%s;time=1235615572.52115;st=coding;_format=Text',...                
    url=sprintf('http://www.ensembl.org/%s/Transcript/Export?db=core;output=fasta;strand=feature;t=%s;st=coding;_format=Text',...                
                sptxt,transid);
	
   %url=sprintf('format=fasta&l=%s%%3A%d-%d&action=export&_format=Text&output=txt&submit=Continue+%%3E%%3E',...
    %url=strcat(base, url);
        
	[s0,status0]=urlread(url);
	if ~(status0>0),
		disp('Unable to download data.');
		return;
	end
	[s] = strread(s0,'%s','delimiter','\n');
	s(1)=[];
	s=char(s);
    
    strs='';
    for k=1:size(s,1)
        strs=[strs,s(k,:)];
    end
    strs=deblank(strs);
