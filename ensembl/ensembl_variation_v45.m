function [rsid,issyn]=ensembl_variation_v45(geneid)

geneid='ENSG00000170095';
rsid={}; issyn=[];
%transid='ENST00000005178';
url=sprintf('http://jun2007.archive.ensembl.org/Homo_sapiens/genesnpview?db=core;gene=%s',geneid);
        
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


mt = find(cellfun(@(y) strcmp(y,'<th class="bottom-border right-border">Validation</th>'),s));
n=length(mt);

if (n>1)
	N=mt(2);
	else
	N=length(s)-1;
end

c=1;
for (k=mt(1)+1:N),
	x=s{k};
        [mat,idx]=regexp(x,'>rs\d+');
	if ~(isempty(mat))
		currentrsid=x(mat+1:idx);
	end        
        if ~isempty(strfind(x,'>SYNONYMOUS_CODING'))
            issyn=[issyn,1];
	    rsid{c}=currentrsid;
	    c=c+1;
        elseif ~isempty(strfind(x,'>NON_SYNONYMOUS_CODING'))
            issyn=[issyn,0];
	    rsid{c}=currentrsid;
	    c=c+1;
        end

end

