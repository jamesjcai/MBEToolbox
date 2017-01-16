function [gid]=genenameapproved(genename)
gid='';
%genename='MTP';
%gid='MTTP';


urlFetch=sprintf('http://www.genenames.org/cgi-bin/hgnc_search.pl?field=all_text&anchor=contains&match=%s&symbol_search=Search&number=5000&format=text&sortby=symbol',...
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
for linum=numLines1+1:numLines2-1
    theline=regexprep(fetchResults{linum},',','');    
    x1=strread(theline,'%s');
    if find(strcmp(genename,x1))>0
        gid=x1{2};
        return;
    end
 end
end
