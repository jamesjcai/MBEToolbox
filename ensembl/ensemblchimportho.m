function [g2]=ensemblchimportho(g1,refspecies)

g2='';
if nargin<3    
    refspecies='Homo_sapiens';
    % refspecies='Gasterosteus_aculeatus';
end

%g1='ENSG00000139618'; g2='ENSPTRG00000005766';
%refspecies='Gasterosteus_aculeatus'; g1='ENSGACG00000016426'; g2='ENSORLG00000000305';

urlFetch=sprintf('http://jun2007.archive.ensembl.org/%s/geneview?gene=%s',...
     refspecies,g1);
try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    rethrow(lasterror);
end



%http://jun2007.archive.ensembl.org/Homo_sapiens/alignview?class=Homology;g
%ene=ENSG00000128573;g1=ENSPTRG00000019608

fetchResults = strread(pagecontent,'%s','delimiter','\n','whitespace','');

for k=1:length(fetchResults)
    x=fetchResults{k};
    [mat,idx] = regexp(x,'g1=ENSPTRG\d+');
if ~isempty(mat)
    g2=x(mat(1)+3:idx(1));
end
end
