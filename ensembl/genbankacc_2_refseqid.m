function [n]=genbankacc_2_refseqid(geneid)

%geneid='U36764';
n='';
try
    urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/nuccore/%s',geneid);
    pagecontent=urlread(urlFetch);

    fetchResults = textscan(pagecontent,'%s','delimiter','\n','whitespace','','bufsize',91256);
    fetchResults = fetchResults{1};
    numLines = strfind(fetchResults,'See reference mRNA sequence for the');
    numLines=find(~cellfun(@isempty,numLines));
    
    if ~(isempty(numLines))
        theline=fetchResults{numLines(1),:};
        theline=strtrim(theline);
        theline=regexprep(theline,'<.*?>','\t');
        x=strfind(theline,'(NM_');
        y=strfind(theline,')');
        n=theline(x+1:y-1);
    end
catch ME
    return;
end