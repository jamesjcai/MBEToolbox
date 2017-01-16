function [model] = modelgy94bgc(omega,kappa,fws,freq)
%MODELGY94BGC - modified GY94 codon model with BGC
%
% Syntax: [model] = modelgy94bgc(omega,kappa,freq)
%
% Inputs:
%    omega   - ratio of nonsyn and syn sub. rates
%    kappa   - $\kappa = \alpha /\beta$
%    freq    - equilibrium frequencies of codons
%
% Outputs:
%    model.R      - Rate matrix
%    morel.freq   - Equilibrium frequency parameters, 1x4
%
% See also:
%REF: Goldman N and Yang Z (1994) MBE

% Molecular Biology & Evolution Toolbox, (C) 2009
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



MAXIMUM_OMEGA = 100;
MAXIMUM_KAPPA = 100;
MAXIMUM_FWS = 1;

MINIMUM_OMEGA = 0.000000;
MINIMUM_KAPPA = 0.000001;
MINIMUM_FWS = 0;

DEFAULT_KAPPA = 2;
DEFAULT_OMEGA = 1;
DEFAULT_FWS = 0.6;

freq=(1/61)*ones(1,61);

if (nargin<1)
	freq=(1/61)*ones(1,61);
	kappa=DEFAULT_KAPPA;
	omega=DEFAULT_OMEGA;
	fws=DEFAULT_FWS;
end


icode=1;
[TABLE,CODON] = codontable;
stops=find(TABLE(icode,:)=='*');
TABLE=TABLE(icode,:);
CODON(stops,:)=[];
TABLE(stops)=[];
TS = [0, 0, 1, 0;
      0, 0, 0, 1;
      1, 0, 0, 0;
      0, 1, 0, 0];

%TV = [0, 1, 0, 1;
%      1, 0, 1, 0;
%      0, 1, 0, 1;
%      1, 0, 1, 0];


TWS = [0, 0, 0, 0;
       1, 0, 0, 1;
       1, 0, 0, 1;
       0, 0, 0, 0];

TSW = [0, 1, 1, 0;
       0, 0, 0, 0;
       0, 0, 0, 0;
       0, 1, 1, 0];


R=zeros(61);
for (i=1:61),
for (j=i:61),
	if (i~=j)
	X=CODON(i,:); Y=CODON(j,:);
	ndiff=sum(X~=Y);
	if (ndiff==1)
		pos=find(X~=Y);
		if(TS(X(1,pos),Y(1,pos))==1)               % transition
			if (TABLE(i)~=TABLE(j))            % transition - nonsynony
			   R(i,j)=kappa*omega;
			else                               % transition - synony
			   R(i,j)=kappa;
			end
		else					   % transversion
			if (TABLE(i)~=TABLE(j))            % transversion - nonsynony
			   R(i,j)=omega;
			else                               % transversion - synony
			   R(i,j)=1;
			end
		end

		if(TWS(X(1,pos),Y(1,pos))==1)               % transition

		elseif(TSW(X(1,pos),Y(1,pos))==1)               % transition

		end

	end
	end
end
end
R=R+R';

model.name='gy94';
model.R=R;
model.freq=freq;