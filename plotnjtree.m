function [anc,brnlen] = plotnjtree(D,showbrnlen,labels,prtnodes,fontsize)
%PLOTNJTREE - Performs neighbor joining (NJ) on distance data and plot NJ tree
%
% Syntax: [anc,brnlen] = plot_njtree(D,showbrnlen,labels,prtnodes,fontsize)
%
% Inputs:
%    D      - Distance matrix
%
% See also: PLOTUPGMA

% RE Strauss, 5/27/96
%   9/7/99 - miscellaneous changes for Matlab v5.
%   9/24/01 - check diagonal elements against eps rather than zero.


% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $




isDebug=0;
if (nargin<2), showbrnlen = 0; end;
if (nargin<3), labels = []; end;
if (nargin<4), prtnodes = []; end;
if (nargin<5), fontsize = []; end;


[n,m]=size(D);
if n~=m
   error('Matrix is not square.')
elseif n<3
   error('Must have at least 3 taxa.')
else

   N=n;            % save original size


if (isDebug)
   names=labels; %   and names of taxa
     if isempty(names)
      for i=1:N
         names=[names, {['S',num2str(i)]} ];
      end
   end
   if sum(size(names)~=[1,N])
      error('Incorrect number of taxon names.')
   end
end


   D=triu(D,1); %wipe out diagonal and lower left triangle
   D=D+D';       %    and copy upper right triangle there


  taxon = 1:N;                            % Initial terminal-taxon labels
  anc = zeros(1,2*N-1);                   % Ancestor function
  brnlen = anc;                            % Branch lengths

   while n>1 % iterative reduction of taxa to 3

      %step 1
      R=sum(D);% compute total distances from one taxon to all others
      M=(n-2)*D-ones(n,1)*R -R'*ones(1,n);% compute M array
      M=M + diag(NaN*ones(n,1));         %      set diagonals to NaN
      [colmins, rowinds]=min(M);
      [minofM, colind] =min(colmins);
      rowind=rowinds(colind); % min entry of M is at (rowind,colind)
      % so taxa rowind and colind are to be joined
    if (isDebug)
          newvertex=['V',num2str(2*N-n+1)];% new internal vertex name
          disp(' ')
          disp(['Taxa ',names{colind},' and ',names{rowind},...
                ' are joined at ',newvertex])
    end

    anc(taxon(colind)) = 2*N-n+1;                   % Attach taxa to common node
    anc(taxon(rowind)) = 2*N-n+1;

      % step 2
      if (n==2)
	dVr=D(rowind,colind)/2;
	dVc=D(rowind,colind)/2;
      else
	dVr=D(rowind,colind)/2+ (R(rowind)-R(colind))/(2*n-4);
	dVc=D(rowind,colind)-dVr; %compute distances from joined taxa to new vertex

	if (isDebug)
	      disp(['Edge from ',names{colind},' to ',newvertex,...
		    ' has length ',num2str(dVc)])
	      disp(['Edge from ',names{rowind},' to ',newvertex,...
		    ' has length ',num2str(dVr)])
	      disp('Hit enter to continue.')
	end


	      if (dVr<0)                % If negative branch length,
		dVc = dVc+dVr;    % Transfer to adjacent branch
		if (dVc<0)                % If adjacent branch negative,
		  dVc = 0;                %   set equal to zero (lengthing the tree)
		end;
		dVr = 0;                   % Set original negative brnlen equal to zero
	      elseif (dVc<0)
		dVr = dVr+dVc;
		if (dVr<0)
		  dVr = 0;
		end;
		dVc = 0;
	      end;
      end
      	    brnlen(taxon(colind)) = dVc;                   % Attach taxa to common node
	    brnlen(taxon(rowind)) = dVr;

      % step 3
    if (isDebug)
          names([colind,rowind])=[];  %remove names of joined taxa
          names=[names,{newvertex}]; %    and add new internal vertex name
    end

      newD=D;
      newD([colind,rowind],:)=[]; %eliminate rows with entries from joined taxa
      newD=[newD, (newD(:,colind)+newD(:,rowind)-D(rowind,colind))/2]; % compute distances to new vertex
      newD(:,[colind,rowind])=[]; % eliminate columns with joined taxa
      n=n-1;
      D=[newD;newD(:,n)' 0];

      taxon([colind,rowind]) = [];
      taxon = [taxon 2*N-n];

   end

   if (showbrnlen)
   %anc
	plot_tree(anc,brnlen,labels,prtnodes,fontsize)
   else
	plot_tree(anc,[],labels,prtnodes,fontsize);
   end

end