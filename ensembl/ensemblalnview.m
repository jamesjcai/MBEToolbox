function [aln]=ensemblalnview(g1,g2,refspecies)

if nargin<3
    refspecies='Homo_sapiens';
    % refspecies='Gasterosteus_aculeatus';
end

%g1='ENSG00000139618'; g2='ENSPTRG00000005766';
%refspecies='Gasterosteus_aculeatus'; g1='ENSGACG00000016426'; g2='ENSORLG00000000305';

urlFetch=sprintf('http://jun2007.archive.ensembl.org/%s/alignview?class=Homology&gene=%s&g1=%s&seq=DNA&format=phylip&submit=Update',...
    refspecies,g1,g2);
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end


fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);

numLines1 = strmatch('</table><pre>',fetchResults);
numLines2 = strmatch('</pre>',fetchResults);

if ~(numLines1>0 && numLines2>0)
    aln=[];
    return;
else

    
linestr=fetchResults{numLines1};
idx=strfind(linestr,'>');
linestr=linestr(idx(2)+1:end);

filename=tempname;
fid=fopen(filename,'w');

    fprintf(fid,'%s\n',linestr);
    
for k=numLines1+1:numLines2-1
    fprintf(fid,'%s\n',fetchResults{k});
end
fclose(fid);

aln = readphylip_i(filename,2,1);
aln.seqnames={g1,g2};
    
    
end

