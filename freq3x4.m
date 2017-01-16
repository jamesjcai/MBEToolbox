function [F] = freq3x4(S)
%FREQ3X4 -

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (isstruct(S)), S=S.seq; end

nn=[];
[n,m]=size(S);
	for (pos=1:3),
		picker=pos:3:m;
    		P=S(:,picker);
		N=[sum(P==1,2),sum(P==2,2),sum(P==3,2),sum(P==4,2)];
		nn=cat(2,nn,N);
	end
m=m/3;
F=nn./m;

if (nargout<1)
      disp('Codon position x base (3x4) table for each sequence.')
      disp(' ')
	for (k=1:n),
	disp(sprintf('#%d',k))
		for (j=0:2),
			disp(sprintf('position  %d:\tA:%.5f\tC:%.5f\tG:%.5f\tT:%.5f',...
			j+1,F(k,1+j*4),F(k,2+j*4),F(k,3+j*4),F(k,4+j*4)))
		end
	disp(' ')
	end

      disp('Codon position x base (3x4) table, overall')
      disp(' ')
	X=sum(nn,1)./(m*n);
	k=1;
	for (j=0:2),
		disp(sprintf('position  %d:\tA:%.5f\tC:%.5f\tG:%.5f\tT:%.5f',...
		j+1,X(k,1+j*4),X(k,2+j*4),X(k,3+j*4),X(k,4+j*4)))
	end
	disp(' ')
end