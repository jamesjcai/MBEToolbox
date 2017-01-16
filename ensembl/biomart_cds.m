function [seq,s]=biomart_cds(transid,species)
%[seq]=biomart_cds(transid,species)
% species = {'human','chimp','macaca'};

if nargin<2, species='human'; end
%transid='ENSPTRT00000000013';
%transid='ENSMMUT00000048360';
%transid='ENST00000003302';

switch species
    case 'chimp'
        url=sprintf('http://www.biomart.org/biomart/martview?VIRTUALSCHEMANAME=default&ATTRIBUTES=ptroglodytes_gene_ensembl.default.sequences.ensembl_transcript_id|ptroglodytes_gene_ensembl.default.sequences.coding&FILTERS=ptroglodytes_gene_ensembl.default.filters.ensembl_transcript_id."%s"&VISIBLEPANEL=resultspanel',...
            transid);
    case 'human'        
        url=sprintf('http://www.biomart.org/biomart/martview?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.coding&FILTERS=hsapiens_gene_ensembl.default.filters.ensembl_transcript_id."%s"&VISIBLEPANEL=resultspanel',...
            transid);
    case 'macaca'
        url=sprintf('http://www.biomart.org/biomart/martview?VIRTUALSCHEMANAME=default&ATTRIBUTES=mmulatta_gene_ensembl.default.sequences.ensembl_transcript_id|mmulatta_gene_ensembl.default.sequences.coding&FILTERS=mmulatta_gene_ensembl.default.filters.ensembl_transcript_id."%s"&VISIBLEPANEL=resultspanel',...
            transid);        
end

try
    [pagecontent,status]=urlread(url);
catch
	disp('Unable to download data.');
	return;
end

yes=false;
fetchResults = strread(pagecontent,'%s','delimiter','\n','whitespace','');
%for k=18968:length(fetchRestults)
c=1;
s='';
for k=18968:length(fetchResults)
    theline=fetchResults{k};
    [mat1]=regexp(theline,'pre class="mart_');
    if mat1>0, yes=true; end
    [mat2]=regexp(theline,'/pre');
    if mat2>0, yes=false; end
    if yes
        if c>1, s=strcat(s,theline); end
        c=c+1;
    end
end
seq=i_encode_n(s);