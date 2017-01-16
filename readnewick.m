function [tree,n,branches,names] = readnewick(file)
%READNEWICK - Reads tree from standart file in Newick format
%
% USAGE: [tree,n,branches,names] = readnewick(file)
%
%
%Adopted from: PHYLLAB toolbox v1.1
%
% function [tree,n,branches,names] = read_newick(tstr)
% Reads tree from standart file in Newick format
%
% Description of this format can be found at
% http://evolution.genetics.washington.edu/phylip/newicktree.html
%
% tree -- a string with the tree
% n -- the number of OTUs
% branches -- branch length
% names -- OTU names
% If no branch length provided, thay are set to ones.
%
%
% Adopted from: PHYLLAB toolbox v1.1

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if nargin < 1
    [file, pathname] = uigetfile( ...
       {'*.tree;*.tre;*.dnd', 'Tree Format Files (*.tree, *.tre, *.dnd)';
        '*.*',  'All Files (*.*)'}, ...
        'Select Phylogenetic Tree File');
	if ~(file), return; end
	file=[pathname,file];
end

% check input is char
% in a future version we may accept also cells
if ~ischar(file)
    error('mbetoolbox:InvalidInput','Input must be a character array')
end

if ~(exist(file,'file') || exist(fullfile(cd,file),'file')),
%  is a valid filename ?
    error('mbetoolbox:InvalidInput','Input must be a valid file')
end





% opening file
ff = fopen(file,'r');
if ff == -1
   n=-1;	tree = [];
   disp('File open error in read_newick.');
   return;
end;

% The tree can be recorded in multiple lines. Removing LF, CR, etc...
[tr, k] = fscanf(ff,'%c');
n=0;
i=1;
while i<k
    if tr(i) <= ' '
       tr(i) =[];
       k=k-1;
    else
       i=i+1;
    end
    if tr(i) == ')'  n=n+1; end
end
n=n+1;

% Now replacing all symbolic names to numbers!
% Even if the name is already number it will be replaced
posO = 1; % position of last '(' - to determin the start of the nest
flagO=0; % the status of previous symbol: 1 - if punctuation, 2 - other
nm = char(n,n);
cur_tip=1;

[l,k] = size(tr);
i=1;
while i<k
	if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':' ) & flagO == 2
      if cur_tip > 1 nm=strvcat(nm,tr(posO+1:i-1));
      else nm=tr(posO+1:i-1);
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
b=zeros(2*n-1,1);
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

      % split string by the ',' symbol
      [token,rem] = strtok(sbst(2:end-1),',');
		rem=rem(2:end);

      % is there a ':' symbol before the ','?
      [t,r] = strtok(token,':');r=r(2:end);
      pr=str2num(t);
      if isempty(r)
         b(pr)=1;
		else
         b(pr)=str2num(r);
      end;

      [t,r] = strtok(rem,':');r=r(2:end);
      vt=str2num(t);
      if isempty(r)
         b(vt)=1;
		else
         b(vt)=str2num(r);
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
   if f t(j)=tree(i); j=j+1; end;
end;

tree=t; branches=b; names=nm;
fclose(ff);