function [n]=pubmed_article_num(geneid)

n=0;
try
    urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/pubmed?term=%s',geneid);
    pagecontent=urlread(urlFetch);


    fetchResults = textscan(pagecontent,'%s','delimiter','\n','whitespace','','bufsize',91256);
    fetchResults = fetchResults{1};
    numLines = strfind(fetchResults,sprintf('articles about %s',geneid));
    numLines=find(~cellfun(@isempty,numLines));
    if isempty(numLines)
        numLines = strfind(fetchResults,sprintf('article about %s',geneid));
        numLines=find(~cellfun(@isempty,numLines));        
    end
    
    if ~(isempty(numLines))
        theline=fetchResults{numLines-1,:};
        theline=strtrim(theline);
        [~,n]=strread(theline,'%s%d');
        % n = regexprep(theline,'<.*?>','\t');
    end
catch ME
    return;
end