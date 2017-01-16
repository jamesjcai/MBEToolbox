function pairs = string2pairs(tree,n)
% function pairs = string2pairs(tree)
% Converts string record of tree to a pair one.

[k,l] = size(tree);
tr = tree;
pairs = zeros(n*2,1); % format from the VOSTORG package

anc = n+1;
posO = 1;
flagO=0;

% parsing the input
i=1;
while  i<=l 
    if tr(i) == '(' 
        flagO = 1;posO = i;
    end;
    if tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ','
        flagO = 2;
    end
    if abs(tr(i))<abs(' ')
        break; 
    end;
    if tr(i) == ')' & flagO == 2
		% now we process the nest
      if anc==2*n break; end;
		sbst = tr(posO:i);
		[pr, np] = sscanf(sbst,'(%d,%d)');
      pairs(pr(1)) = anc;         
      pairs(pr(2)) = anc;                       
      % now we rewrite the tree string
      tr = sprintf('%s%d',tr(1:posO-1),anc,tr(i+1:l));
      anc=anc+1;
      i=0;
      [k,l]=size(tr);
   end;
	i=i+1;     
end;
pairs(anc-1)=anc;   
