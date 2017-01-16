function [ensmid,geneid,genename2]=genenamesearchX(genename)

ensmid='';
geneid='';
genename2='';

[ensmid,geneid]=genenamesearch3(genename);

if isempty(ensmid)
    genename2=genenameapproved(genename);
    if isempty(genename2)
        genename2=genenameapproved(upper(genename));
        if isempty(genename2)
            return;
        else
            [ensmid,geneid]=genenamesearch3(genename2);                    
        end        
    else
            [ensmid,geneid]=genenamesearch3(genename2);        
        return;
    end
else
    return;
end

