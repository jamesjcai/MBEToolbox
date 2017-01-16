function [n]=pubmed_article_num2(geneid)

n=0;
try
    urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/gene?term=(%s%%5BGene%%20Name%%5D)%%20AND%%20Homo%%20sapiens%%5BOrganism%5D',geneid);
    pagecontent=urlread(urlFetch);

    fetchResults = textscan(pagecontent,'%s','delimiter','\n','whitespace','','bufsize',91256);
    fetchResults = fetchResults{1};
    numLines = strfind(fetchResults,'See all');
    numLines=find(~cellfun(@isempty,numLines));
    
    if ~(isempty(numLines))
        theline=fetchResults{numLines,:};
        theline=strtrim(theline);
        theline=regexprep(theline,'<.*?>','\t');
        x=strfind(theline,'See all (');
        y=strfind(theline,') citation');
        n=str2double(theline(x+9:y-1));        
    end
catch ME
    return;
end