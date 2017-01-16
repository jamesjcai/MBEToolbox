function [fid]=ensemblg2f(gid,speciesid)

%gid=ENSG00000177102
%e.g., gid='ENSG00000137975', fid=fam50v00000001049

if nargin<2, speciesid=1; end
%e.g., pid=ENSP00000324595, gid=ENSG00000177102

spename={'Homo_sapiens','Pan_troglodytes',...
'Mus_musculus','Rattus_norvegicus','Canis_familiaris'};

%tagname={'ENSG','ENSPTRG',...
%'ENSMUSG','ENSRNOG','ENSCAFG'};

%ptagname={'ENSP','ENSPTRP',...
%'ENSMUSP','ENSRNOP','ENSCAFP'};


fid='';
urlFetch=sprintf('http://jul2008.archive.ensembl.org/%s/geneview?gene=%s',spename{speciesid},gid);

%urlFetch
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end

%      <a href="/Homo_sapiens/familyview?family=fam50v00000001049">fam50v00000001049</a> : CALCIUM ACTIVATED CHLORIDE CHANNEL 

fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
numLines = strfind(fetchResults,'fam50v');
numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);

if ~(isempty(numLines))
    theline=fetchResults(numLines,:);
    theline=theline{1};
   
    %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
    %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
    
    %theline=strrep(theline,'X','23');
    %theline=strrep(theline,'Y','24');
[mat, idx] = regexp(theline,'fam50v\d\d\d\d\d\d\d\d\d\d\d');

    fid=theline(mat:idx);
end
