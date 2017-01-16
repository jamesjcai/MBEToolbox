function c = mbelfcreator(tree, n)
%  c = mbelfcreator(tree, n, number) outputs source of
%                      a MatLab subroutine  
%                      for computing a likelihood value for one
%                      site of n homologous sequences:
%                            *s(1) ..., s(n)*
%                      P is array of transition matrices
%                      times *t(1), ..., t(2n-1)*
%                      freq - a root distribution (IMPORTANT -- ROW vector!)
%                      (note that we don't need to know the values of
%                      s and t to get the commands);
%
%      note that each generated subroutine can have a unique name due to 
%      addition of *number* to the name of the new subroutine
%
% PHYLLAB toolbox v1.0 internal function

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


% step 0: initialize
    ancestor = n+1;     % current ancestral node
    ncommand = 1;       % current command to store
    nodes = zeros(n,1); % current set of children
                        % nnodes is the number of children

    ss = cell(2*n-1);
    for i=1:2:n-1
       ss{i}=cellstr('');
    end
% generates a subroutine for computing likelihood (for one site) for a tree 
% with n terminal branches
    eol = [setstr(13) setstr(10)];
    
    %s = sprintf('LL = lhi(P,n,rou,s,nalpha)');
    
    c = [];
    [y x] = size(tree);
    for i = 1:x
        c1 = tree(i);
        if strcmp(c1,',')==1
           tree(i) = ' '; 
        end
    end
% unpack matrices from the common array and truncate those for terminal brunches    
    for i = 1:n
        c = sprintf('P((%d*nalpha+1):(%d*nalpha),s(%d))',i-1,i,i);
	%c = sprintf('P%d',i);
        ss{i} = cellstr(c); 
    end

    for i = n+1:2*n-2					% James J. Cai		n-3 --> n-2
       c = sprintf('P((%d*nalpha+1):(%d*nalpha),:)',i-1,i); 
       %c = sprintf('P%d',i)
       ss{i} = cellstr(c);
    end
      c = sprintf('eye(nalpha)');
      ss{2*n-1} = cellstr(c);				 % James J. Cai		n-2 --> n-1
% step # ancestor - n

for ancestor = n+1:2*n-1
   
        [l,r,s] = findnest(tree);           % found the next nested set
                                            % multifurcating tree
                  if l == 0 | r == 0
                            ancestor = ancestor - 1;
                            break;
                  end
	          if ancestor > n+1
                    nbran = ancestor-1;
                  end
        [y,x] = size(s);
        s=s(2:x);
        [nodes, nnodes] = sscanf(s,'%d');   % determined the list of OTUs
        
        %disp(nodes);
        %disp(nnodes);
        
        c = '(';
        
        % Pi = Pn1.* ... .*Pnx;
        for i = 1:nnodes
           
          if i < nnodes
              c = [c char(ss{nodes(i)}) '.*'];        
          else
              c = [c char(ss{nodes(i)})];
          end
        end
        
        
        if strcmp(char(ss{ancestor}),'') ~= 1 
           c = ['(' char(ss{ancestor}) '*' c '))' ];
        else
           c = [c ')'];
        end
        
        %disp(c);    
        ss{ancestor}=cellstr(c);
        [y,mm] = size(tree);
        tr1 = ' '; tr2 = ' ';
        if l > 1 & r < mm
             tr1 = tree(1:l-1);
             tr2 = tree(r+1:mm);
        end
        s = num2str(ancestor);
        tree = [tr1 s tr2];
end
  c = ['(freq*' char(ss{ancestor}) ')'];


function [l,r,s1] = findnest(s)
% [l,r,s1] = findnest(s), where s - is a "tree string", 
% finds the position of the leftmost terminal cluster
% l & r are left and right edges of substring
% s1 is the substring found           
%
% PHYLLAB toolbox v1.0 internal function
l = 0; r = 0;
[l,k] = size(s); 
c = zeros(1,k); 
nadr = 1; 
address = zeros(1,k);
% 2, 3 are left and right brackets, respectively
le = '('; ra = ')';
count = 0;
for i = 1:k
    if     strcmp(s(i),le) == 1
           c(nadr) = 2;
           address(nadr)=i;
           nadr = nadr + 1;
           count = count+1;
    elseif strcmp(s(i),ra) == 1
           c(nadr) = 3; 
           address(nadr)=i;
           nadr = nadr + 1;
           count = count+1;
    end 
end 
if count < 2
		r=0; %%%%%%%%%%%%%%%%%%
		s1 = s; %%%%%%%%%%%%%%%
      return;
end
for i = 2:nadr
    if  c(i-1) == 2 & c(i) == 3
        l = address(i-1);
        r = address(i);
        break;
    end
end
address = []; c = [];        
s1 = s(l:r);