function [gid]=genenamesearch2(genename)
gid='';
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
    [gidx]=i_getensg(urlFetch);
    if ~isempty(gidx)
        gid=gidx;
    end
end

%web(urlFetch)

function [gid]=i_getensg(urlFetch)
gid='';

try
    pagecontent=urlread(urlFetch);

fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);

numLines1 = strmatch('<td>ENSG',fetchResults);
if ~(isempty(numLines1))
    theline=fetchResults{numLines1};
    [mat1,mat2] = regexp(theline,'\d+');
    gid=theline(mat1-4:mat2);
end
catch
    %errordlg(lasterr)
    disp(urlFetch)
    
    rethrow(lasterror);
end


