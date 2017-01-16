function [rsid,issyn,alletype]=ensembl_variation(transid)



rsid={}; issyn=[]; alletype={};
%transid='ENST00000005178';
url=sprintf('http://www.ensembl.org/Homo_sapiens/Component/Transcript/Web/ProteinVariations?db=core;t=%s;vdb=variation',transid);
        
	[s0,status0]=urlread(url);
	if ~(status0>0),
		disp('Unable to download data.');
		return;
	end
	[s] = strread(s0,'%s','delimiter','\n');    
    

if length(s)<10, return; end
while isempty(s{1})
    s(1) = [];
end

mt = find(cellfun(@(y) strcmp(y,'</tr>'),s));
n=length(mt)-1;

if n>0
    issyn=ones(1,n);
    for k=1:n
        x=s{mt(k)+3};
        [mat,idx] = regexp(x,'>rs\d+');
        rsid{k}=x(mat+1:idx);
        x=s{mt(k)+4};
        if ~isempty(strfind(x,'Synonymous'))
            issyn(k)=1;
        elseif ~isempty(strfind(x,'Non-synonymous'))
            issyn(k)=0;
        else
            issyn(k)=nan;
        end
        if nargout>2
            x=s{mt(k)+5};        
            [mat,idx] = regexp(x,'>\D+<');
            alletype{k}=x(mat+1:idx-1);     
        end
    end    
end
    