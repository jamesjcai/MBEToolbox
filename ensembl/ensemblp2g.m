function [gid]=ensemblp2g(pid,speciesid)


if nargin<2, speciesid=1; end
%e.g., pid=ENSP00000324595, gid=ENSG00000177102

spename={'Homo_sapiens','Pan_troglodytes',...
'Mus_musculus','Rattus_norvegicus','Canis_familiaris'};

tagname={'ENSG','ENSPTRG',...
'ENSMUSG','ENSRNOG','ENSCAFG'};

%ptagname={'ENSP','ENSPTRP',...
%'ENSMUSP','ENSRNOP','ENSCAFP'};


gid='';
urlFetch=sprintf('http://www.ensembl.org/%s/protview?peptide=%s',spename{speciesid},pid);
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end
    
fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
numLines = strfind(fetchResults,tagname{speciesid});
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);

if ~(isempty(numLines))
    theline=fetchResults(numLines,:);
    theline=theline{1};
   
    %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
    %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
    
    %theline=strrep(theline,'X','23');
    %theline=strrep(theline,'Y','24');
[mat,idx] = regexp(theline,...
    sprintf('%s\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d',tagname{speciesid}));

    gid=theline(mat:idx);
end
