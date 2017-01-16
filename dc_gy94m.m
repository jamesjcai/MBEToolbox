function [dS,dN,dN_dS,lnL,para] = dc_gy94m(aln,a,b)
%DC_GY94M - dS, dN estimation by modified GY94 model
%
% [dS,dN] = dc_gy94m(aln,a,b)
% calculates dS and dN between sequence a and b in aln.
%
%%

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if nargin < 3, error('DC_ML requires at least three input arguments');  end

global noise
noise=1;

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[aln]=rmcodongaps(seq);
s1=seq(a,:); s2=seq(b,:);

%if (nargin<4)
%	kappa=1.6;   % fixed kappa
%	omega=0.8;
%	md=modelgy94m(omega,kappa);
%end

%%
% Guess: if not codonise61ed then do it
%%
if (sum(s1>5)<2 & sum(s2>5)<2),
	s1=codonise61(s1); s2=codonise61(s2);
	% disp('s1 and s2 have been codonise61 now')
end


%select=nargin-2;
select=1;
switch (select)
    case (1)
	disp('optimising everything: t, kappa1, kappa2 and omega')
	[para,lnL]=i_optimtk2o(s1,s2);
	t=para(1);
	kappa1=para(2);
	kappa2=para(3);
	omega=para(4);
	md=modelgy94m(omega,kappa1,kappa2);		% build model from optimised values
    case (2)
    	disp('fixed kappa, optimising t and omega.')
	error('under development!')
    case (3)
    	disp('fixed kappa and omega, optimising t.')
	error('under development!')
	kappa1=1.6;   % fixed kappa
	kappa2=1.7;   % fixed kappa
	omega=0.8;
	md=modelgy94m(omega,kappa1,kappa2);
	[t,lnL] = optimlikelidist(md,s1,s2,0,2);
    otherwise
        error('invalid selection!')
end




%%
% Composes substitution rate matrix, Q
%%
Q=composeQ(md.R,diag(md.freq))./61;
%[V,D] = eig(Q*t);
%P=V*diag(exp(diag(D)))/V;
%%P=expm(Q*t);


%%
% Making a mask matrix, M
%%
icode=1; [TABLE] = codontable;
stops=find(TABLE(icode,:)=='*');
TABLE=TABLE(icode,:);
TABLE(stops)=[];
M=zeros(61);
for (i=1:61),
for (j=i:61),
	if (i~=j)
	if (TABLE(i)==TABLE(j))     % synony changes
	   M(i,j)=1;
	end
	end
end
end
M=M+M';


%%
% Calculate pS and pN, when omega = optimised omega
%%
pS=sum(sum(Q.*M));
pN=1-pS;
pS=pS*t;     % why *t here?
pN=pN*t;


%%
% Calculate pS and pN when omega = 1
%%
omega2=1;
md2=modelgy94m(omega2,kappa1,kappa2);
Q2=composeQ(md2.R,diag(md2.freq))./61;
pS2=sum(sum(Q2.*M));
pN2=1-pS2;
pS2=pS2*3; pN2=pN2*3;

%%
% Calculates dS and dN
%%
dS=pS/pS2;
dN=pN/pN2;

if (nargout>2), dN_dS=dN./dS; end


disp(' ')
disp('NOTE: Here we used an equal codon frequence!')
disp('(i.e., in codeml setting: CodonFreq = 0  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table)')
disp(' ')


function [para,lnL]=i_optimtk2o(s1,s2)
	global noise;
	et=0.5; ek1=1.5; ek2=1.2; eo=0.8;    % initial values for t, kappa and omega
	options = optimset('fminsearch');

	if (noise), options=optimset(options,'display','iter');
	else, options=optimset(options,'display','off');
	end

	[para,f_opt]=fminsearch(@i_likelifuntk2o,[et,ek1,ek2,eo],options,s1,s2);
	lnL=-f_opt;
	if (noise),
		disp(sprintf('lnL = %.5f',lnL))
		disp(sprintf('t=%.5f, kappa=%.5f, omega=%.5f',para(1),para(2),para(3)))
	end

function [lnL] = i_likelifuntk2o(x,s1,s2)
	lnL=inf;
	if (any(x<eps)|any(x>20)), return; end
	t=x(1); kappa1=x(2); kappa2=x(3); omega=x(4);
	if (t<eps|t>5), return; end
	if (kappa1<eps|kappa1>20), return; end
	if (kappa2<eps|kappa2>20), return; end
	if (omega<eps|omega>10), return; end
	md=modelgy94m(omega,kappa1,kappa2);
	[lnL] = -1*likelidist(t,md,s1,s2);

