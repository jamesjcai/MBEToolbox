function [gid,pagecontent]=ensemblf2g(fid,speciesid)

%gid=ENSG00000177102

%e.g., fid='fam50v00000001049', gid=ENSG00000137975

if nargin<2, speciesid=1; end
%e.g., pid=ENSP00000324595, gid=ENSG00000177102

spename={'Homo_sapiens','Pan_troglodytes',...
'Mus_musculus','Rattus_norvegicus','Canis_familiaris','Macaca_mulatta'};

tagname={'ENSG','ENSPTRG',...
'ENSMUSG','ENSRNOG','ENSCAFG','ENSMMUG'};

%ptagname={'ENSP','ENSPTRP',...
%'ENSMUSP','ENSRNOP','ENSCAFP'};


gid={''};
urlFetch=sprintf('http://jul2008.archive.ensembl.org/%s/familyview?family=%s',...
    spename{speciesid},fid);
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end

%   <a href="/Homo_sapiens/familyview?family=fam50v00000001049">fam50v00000001049</a> : CALCIUM ACTIVATED CHLORIDE CHANNEL 
fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);


%{
hstartn=0; hendn=0;
pstartn=0; pendn=0;

numLines = strfind(fetchResults,'Location of Ensembl genes containing family');
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, hstartn=numLines(1); end

numLines = strfind(fetchResults,'Other peptides in Family');
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, hendn=numLines(1); end

numLines = strfind(fetchResults,'Pan troglodytes');
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, pstartn=numLines(1); end

numLines = strfind(fetchResults,'ENSPTRP');
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, pendn=numLines(end); end
%}


[hstartn,hendn]=i_linenumpicker(fetchResults,...
    'Location of Ensembl genes containing family',...
    'Other peptides in Family');
[gid]=i_ensmidpicker(fetchResults,hstartn,hendn,tagname{speciesid});




function [pstartn,pendn]=i_linenumpicker(fetchResults,targettag1,targettag2)
pstartn=0; pendn=0;

numLines = strfind(fetchResults,targettag1);
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, pstartn=numLines(1); end

numLines = strfind(fetchResults,targettag2);
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);
if numel(numLines)>0, pendn=numLines(end); end



function [gid2]=i_ensmidpicker(fetchResults,pstartn,pendn,targettag)

gid2={''};
if (pstartn>0&&pendn>0)    
    c=0;
    thelines=fetchResults(pstartn:pendn,:);
    
    for k=1:length(thelines)
    theline=thelines{k};   
    %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
    %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
    
    %theline=strrep(theline,'X','23');
    %theline=strrep(theline,'Y','24');
    [mat,idx] = regexp(theline,...
        sprintf('%s\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d',targettag));
    
    if ~isempty(theline(mat:idx))
        c=c+1;
        gid2{c}=theline(mat:idx);
    end
    end
    gid2=unique(gid2);
end

