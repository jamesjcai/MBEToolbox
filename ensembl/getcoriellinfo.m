function [txt]=getcoriellinfo(id)

%id='NA12236';
%id='NA18537';

txt='';
try
    urlFetch=sprintf('http://ccr.coriell.org/Sections/Search/Search.aspx?PgId=165&q=%s',id);
    pagecontent=urlread(urlFetch);


    fetchResults = textscan(pagecontent,'%s','delimiter','\n','whitespace','','bufsize',91256);
    fetchResults = fetchResults{1};
    numLines = strfind(fetchResults,'Sample_Detail');
    numLines=find(~cellfun(@isempty,numLines));
    %numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);

    if ~(isempty(numLines))
        theline=fetchResults{numLines,:};
        %theline=theline{1};
        %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
        %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
        %theline=strrep(theline,'X','23');
        %theline=strrep(theline,'Y','24');
        %[mat, idx] = regexp(theline,'fam50v\d\d\d\d\d\d\d\d\d\d\d');
        txt = regexprep(theline,'<.*?>','\t');

    end
catch ME
    txt='';
end