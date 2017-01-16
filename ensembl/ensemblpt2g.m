function [gid]=ensemblpt2g(ptid)

% see also: genenamesearch3

gid='';
[temp,a]=regexp(ptid,'ENSP','match','start');
if isempty(a)    
    [temp,a]=regexp(ptid,'ENST','match','start');
if isempty(a)    
    return;
else
%urlFetch=sprintf('http://www.ensembl.org/Homo_sapiens/transview?transcript=%s',ptid);        
urlFetch=sprintf('http://www.ensembl.org/Homo_sapiens/Transcript/Transcript?t=%s',ptid);
end
else
urlFetch=sprintf('http://www.ensembl.org/Homo_sapiens/protview?peptide=%s',ptid);    
end

%e.g., pid=ENSP00000324595, gid=ENSG00000177102

%ENST00000252900
%http://www.ensembl.org/Homo_sapiens/transview?db=core;transcript=ENST00000252900




try
     pagecontent=urlread(urlFetch);
catch
     %errordlg(lasterr)
     disp(urlFetch)
     rethrow(lasterror);
end

fetchResults =char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
numLines = strfind(fetchResults,'ENSG'); numLines=find(~cellfun(@isempty,numLines));
%numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);

if ~(isempty(numLines))
     theline=fetchResults(numLines,:);
     theline=theline{1};

     %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
     %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only

     %theline=strrep(theline,'X','23');
     %theline=strrep(theline,'Y','24');
[mat, idx] = regexp(theline,'ENSG\d\d\d\d\d\d\d\d\d\d\d');

     gid=theline(mat:idx);
end

