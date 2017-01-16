function [genename]=omimid2genename(omimid)

%see also genename2omimid
%omimid='603904';
genename='';
if ischar(omimid)
urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s&cmd=plain',...
    omimid);
else
    urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%d&cmd=plain',...
    omimid);
end

try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end


fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
    theline=fetchResults{1};
    [mat1,mat2] = regexp(theline,';');
    genename=strtrim(theline(mat1+1:end));