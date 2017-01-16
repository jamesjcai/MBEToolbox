function [seqtxt2]=deambiguity(seqtxt)
%DEAMBIGUITY - de-ambiguity DNA
%
%

%        R = G A (purine)        
%        Y = T C (pyrimidine)    
%        K = G T (keto)    
%        M = A C (amino)
%        S = G C 
%        W = A T 
%
%        B = G T C
%        D = G A T
%        H = A C T
%        V = G C A
%        N = A G C T (any)  

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(seqtxt);

seqtxt2=[];
for k=1:n
    s=seqtxt(k,:);
    s2=i_builds2(s);
    seqtxt2=[seqtxt2;s2];
end



function s2=i_builds2(s)
    s2=[];
    
    ambiChar='RYKMSW';
    acgtChar={'GA','TC','GT','AC','CC','AT'};

%    ambiChar='BDHV';
%    acgtChar={'GTC','GAT','ACT','GCA'};
%    ambiChar='N';
%    acgtChar={'AGCT'};    

    for k=1:length(ambiChar)
        ambic=ambiChar(k);
        acgtc=acgtChar{k};
    
        idx=find(s==ambic);
        if ~isempty(idx)
            for i=1:length(acgtc)-size(s2,1), s2=[s2;s]; end
            fillings=transpose(acgtc);
            for j=1:length(idx), s2(:,idx(j))=fillings; end
        end    
    end
    if isempty(s2)
        s2=s;
    end
    
    
    