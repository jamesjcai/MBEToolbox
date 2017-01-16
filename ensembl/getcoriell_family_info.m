function getcoriell_family_info(id)

%id='NA12236';
%id='1334';

txt='';
outtxt='';
try
    urlFetch=sprintf('http://ccr.coriell.org/sections/BrowseCatalog/FamilyTypeSubDetail.aspx?PgId=402&fam=%s&coll=GM',id);
    pagecontent=urlread(urlFetch);


    fetchResults = textscan(pagecontent,'%s','delimiter','\n','whitespace','','bufsize',91256);
    fetchResults = fetchResults{1};
    %numLines = strfind(fetchResults,'Catalog ID');
    numLines = strfind(fetchResults,'Sample_Detail.aspx');
    numLines=find(~cellfun(@isempty,numLines));    
    %numLines = strmatch('     <li class=""bullet""> <a href=""/Homo_sapiens/geneview?db=core;gene=',fetchResults);

    if ~(isempty(numLines))
        for ss=1:length(numLines)
        theline=fetchResults{numLines(ss),:};
        %theline=theline{1};
        %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
        %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
        %theline=strrep(theline,'X','23');
        %theline=strrep(theline,'Y','24');
        %[mat, idx] = regexp(theline,'fam50v\d\d\d\d\d\d\d\d\d\d\d');
        txt = regexprep(theline,'<.*?>','\t');
        txt=strtrim(txt);
        txt = regexprep(txt,'\t\t','\t');
        txt = regexprep(txt,'\t\t','\t');
        txt = regexprep(txt,'\t\t','\t');
        %outtxt=[outtxt,'\n',txt];
        fprintf('%s\n',txt);
        end
    end
catch ME
    outtxt='';
end

