function [topo,notu,brchlen,names] = parsetree(tree)
%PARSETREE - Parses string of tree in Newick format
%
% Syntax: [topo,notu,brchlen,names] = parsetree(tree)
%
% [tree,notu,brchlen,names] = parsetree(tree)
% Reads tree from standart file in Newick format
%
% Description of this format can be found at
% http://evolution.genetics.washington.edu/phylip/newicktree.html
%
% tree -- a string with the tree
% notu -- the number of OTUs (n_otu)
% brchlen -- branch length
% names -- OTU names
% If no branch length provided, thay are set to ones.

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


tr=removeblanks(tree);
k=length(tr);

notu=0;
i=1;
while i<k
    if tr(i) <= ' '
       tr(i) =[];
       k=k-1;
    else
       i=i+1;
    end
    if tr(i) == ')'  notu=notu+1; end
end
notu=notu+1;

% Now replacing all symbolic names to numbers!
% Even if the name is already number it will be replaced
posO = 1; % position of last '(' - to determin the start of the nest
flagO=0; % the status of previous symbol: 1 - if punctuation, 2 - other
names = char(notu,notu);
cur_tip=1;

[l,k] = size(tr);
i=1;
while i<k
	if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':' ) & flagO == 2
      if cur_tip > 1 names=strvcat(names,tr(posO+1:i-1));
      else names=tr(posO+1:i-1);
		end;
      tr = sprintf('%s%d%s',tr(1:posO),cur_tip,tr(i:end));
      [j,i]=size(sprintf('%s%d',tr(1:posO),cur_tip));
      i=i+1;
      [j,k]=size(tr);
      cur_tip=cur_tip+1;
   end;
   if (tr(i) == '(' | tr(i) == ')' | tr(i) == ',')
        flagO = 1;posO = i;
   end;
   if tr(i) == ':'
      flagO = 3;
   end;
   if (tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ',' & tr(i) ~= ';' & tr(i) ~= ':' & flagO ~= 3)
      flagO = 2;
   end;
	i=i+1;
end;

% now we will deal with branch length

[k,l] = size(tr);
tree=tr; % save the tree
brchlen=zeros(2*notu-1,1);
anc = notu+1;
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
      if anc==2*notu break; end;
      sbst = tr(posO:i);

      % split string by the ',' symbol
      [token,rem] = strtok(sbst(2:end-1),',');
		rem=rem(2:end);

      % is there a ':' symbol before the ','?
      [topo,r] = strtok(token,':');r=r(2:end);
      pr=str2num(topo);
      if isempty(r)
         brchlen(pr)=1;
		else
         brchlen(pr)=str2num(r);
      end;

      [topo,r] = strtok(rem,':'); r=r(2:end);
      vt=str2num(topo);
      if isempty(r),
         brchlen(vt)=1.0;
      else
         brchlen(vt)=str2num(r);
      end;

      % now we rewrite the tree string
      tr = sprintf('%s%d',tr(1:posO-1),anc,tr(i+1:l));
      anc=anc+1;
      i=0;
      [k,l]=size(tr);
   end;
	i=i+1;
end;

% remove branch length from the tree

[k,l] = size(tree);
f = 1;
j=1;
for i=1:l
   if tree(i) == ':' f = 0; end;
   if tree(i) == ')' | tree(i) == '(' | tree(i) == ','  f = 1; end;
   if f topo(j)=tree(i); j=j+1; end;
end;

if (nargout<1),
	numbrch = 2*notu-2;
	if ((length(brchlen)-1==numbrch)&&(length(names)==notu)),
	 disp('Rooted tree.')
	else
	 disp('Unrooted tree.')
	end
end
