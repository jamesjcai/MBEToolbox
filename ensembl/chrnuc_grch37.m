function [raws,s] = chrnuc_grch37(chrid,startn,endn,species)
%CHRNUC_GRCH37 - returns nucleotides of human geneome
%USAGE: chrnuc_grch37(chrid,startn,endn)

% Population Genetics & Evolution Toolbox, (C) 2010
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2010

if (nargin<3), endn=startn; end
if (nargin<4), species='human'; end

switch lower(species)
    case 'human'
        sptxt='Homo_sapiens';
    case 'mouse'        
        sptxt='Mus_musculus';
    otherwise
        error('Species is not supported.')
end


if (ischar(chrid)), 
	chrid=upper(chrid);
	if (strcmp(chrid,'23')), chrid='X'; end
	if (strcmp(chrid,'24')), chrid='Y'; end
elseif(chrid==23), 
	chrid='X';
elseif(chrid==24), 
	chrid='Y';
else
	chrid=int2str(chrid); 
end



[bp] = chrlen(chrid);
if ~(bp>0), error('Need valid chromosome name, such as ''2'', ''X'''); end
if (startn>endn), error('startn should < endn'); end
%if (endn>bp), error('Out of range.'); end

[s,raws]=i_getnuc(chrid,startn,endn,sptxt);



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [strs,s] = i_getnuc(chorom, startn, endn,sptxt)

    %output=fasta;r=6:133017695-133017699;strand=1;genomic=unmasked;_format=Text
%	base=sprintf('http://www.ensembl.org/%s/Location/Export?',sptxt);	
%	base=sprintf('http://aug2006.archive.ensembl.org/%s/exportview?',sptxt);	    
%	url=sprintf('output=fasta&r=%s%%3A%d-%d&strand=1;genomic=unmasked;_format=Text',...
%                chorom, startn, endn);           

	base=sprintf('http://useast.ensembl.org/%s/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;strand=feature;coding=yes;cdna=yes;peptide=yes;utr3=yes;exon=yes;intron=yes;genomic=unmasked;utr5=yes;_format=Text;',sptxt);
    url=sprintf('r=%s%%3A%d-%d',chorom, startn, endn);
    url=strcat(base, url);
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
