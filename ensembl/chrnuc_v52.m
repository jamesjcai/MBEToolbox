function [raws,s] = chrnuc_v52(chrid,startn,endn,species)
%CHRNUC_V52 - returns nucleotides of reference geneome
%USAGE: chrnuc_v52(chrid,startn,endn)

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007

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

	base=sprintf('http://www.ensembl.org/%s/Location/Export?',sptxt);	
%	base=sprintf('http://aug2006.archive.ensembl.org/%s/exportview?',sptxt);	    
    url=sprintf('output=fasta;r=%s%%3A%d-%d;strand=1;st=genomic;_format=Text',...
                chorom, startn, endn);            
	%url=sprintf('format=fasta&l=%s%%3A%d-%d&action=export&_format=Text&output=txt&submit=Continue+%%3E%%3E',...
    url=strcat(base, url);
    url
        
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
