function [model] = modelgy94m(omega,kappa1,kappa2,freq)
%MODELGY94M - Modified Goldman & Yang (1994) codon model
%
% Syntax: [model] = modelgy94m(omega,kappa1,kappa2,freq)
%
% Inputs:
%    kappa1   -
%
% Outputs:
%    model.R      - Rate matrix
%    morel.freq   - Equilibrium frequency parameters, 1x4
%
% See also:
%REF: Goldman N and Yang Z (1994) MBE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


MAXIMUM_OMEGA = 100;
MAXIMUM_KAPPA1 = 100;
MAXIMUM_KAPPA2 = 100;
MINIMUM_OMEGA = 0.000000;
MINIMUM_KAPPA = 0.000001;

DEFAULT_KAPPA1 = 1.2;
DEFAULT_KAPPA2 = 1.5;
DEFAULT_OMEGA = 1;

freq=(1/61)*ones(1,61);

if (nargin<1)
	freq=(1/61)*ones(1,61);
	kappa1=DEFAULT_KAPPA1;
	kappa2=DEFAULT_KAPPA2;
	omega=DEFAULT_OMEGA;
end


icode=1;
[TABLE,CODON] = codontable;
stops=find(TABLE(icode,:)=='*');
TABLE=TABLE(icode,:);
CODON(stops,:)=[];
TABLE(stops)=[];

TS1= [0, 0, 0, 0;
      0, 0, 0, 1;
      0, 0, 0, 0;
      0, 1, 0, 0];

TS2= [0, 0, 1, 0;
      0, 0, 0, 0;
      1, 0, 0, 0;
      0, 0, 0, 0];

TV = [0, 1, 0, 1;
      1, 0, 1, 0;
      0, 1, 0, 1;
      1, 0, 1, 0];

R=zeros(61);
for (i=1:61),
for (j=i:61),
	if (i~=j)
	X=CODON(i,:); Y=CODON(j,:);
	ndiff=sum(X~=Y);
	if (ndiff==1)
		pos=find(X~=Y);
		if(TS1(X(1,pos),Y(1,pos))==1)               % transition
			if (TABLE(i)~=TABLE(j))            % transition - nonsynony
			   R(i,j)=kappa1*omega;
			else                               % transition - synony
			   R(i,j)=kappa1;
			end
		elseif(TS2(X(1,pos),Y(1,pos))==1)               % transition

			if (TABLE(i)~=TABLE(j))            % transition - nonsynony
			   R(i,j)=kappa2*omega;
			else                               % transition - synony
			   R(i,j)=kappa2;
			end
		else					   % transversion
			if (TABLE(i)~=TABLE(j))            % transversion - nonsynony
			   R(i,j)=omega;
			else                               % transversion - synony
			   R(i,j)=1;
			end
		end
	end
	end
end
end
R=R+R';

model.name='gy94m';
model.R=R;
model.freq=freq;