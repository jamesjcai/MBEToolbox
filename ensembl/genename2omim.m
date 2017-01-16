function [omimid]=genename2omim(genename)

% see also: ensemblp2g

omimid='';

%genename='PDHA1';
%genename='MTM1';

urlFetch=sprintf('http://www.genenames.org/cgi-bin/hgnc_search.pl?field=symbol&anchor=equals&match=%s&symbol_search=Search&number=50&format=text&sortby=symbol',...
    genename);

%urlFetch=sprintf('http://www.genenames.org/cgi-bin/hgnc_search.pl?sortby=symbol&number=100&format=html&field=all_text&andnot=AND&anchor=contains&sortby2=symbol&search_type=simple&match=%s',...
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end

fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);
numLines1 = strmatch('</head><body><pre>',fetchResults);
numLines2 = strmatch('</pre></body></html>',fetchResults);


if ~(isempty(numLines1))
    if (numLines2-numLines1==2)
        %disp('OK');
    end
    
    theline=fetchResults{numLines1+1};
    [mat1,mat2] = regexp(theline,'^\d+');
    
    urlFetch=sprintf('http://www.genenames.org/data/hgnc_data.php?hgnc_id=%s',theline(mat1:mat2));
    [ensmidx]=i_getensg(urlFetch);
    if ~isempty(ensmidx)
        omimid=ensmidx;
    end

end

%web(urlFetch)

function [ensmid]=i_getensg(urlFetch)
ensmid='';


try
    pagecontent=urlread(urlFetch);

fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);

numLines1 = strmatch('<td colspan="4"><a href="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=',fetchResults);
if ~(isempty(numLines1))
    theline=fetchResults{numLines1};
    [mat1,mat2] = regexp(theline,'\d+');
    ensmid=theline(mat1(2):mat2(2));
end

catch
    %errordlg(lasterr)
    disp(urlFetch)
    
    rethrow(lasterror);
end


