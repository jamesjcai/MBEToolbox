function [dS,dN,dN_dS,lnL,value] = dc_gy94(aln,a,b)
%DC_GY - dS, dN estimation by codeml method
%
% [dS,dN] = dc_gy94(aln,a,b)
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



if nargin < 3, error('DC_GY94 requires at least three input arguments');  end

global noise
noise=1;

if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[seq]=rmcodongaps(seq);
s1=seq(a,:); s2=seq(b,:);

%if (nargin<4)
%	kappa=1.6;   % fixed kappa
%	omega=0.8;
%	md=modelgy94(omega,kappa);
%end

%%
% Guess: if not codonise61ed then do it
%%
if (sum(s1>5)<2 && sum(s2>5)<2),
	s1=codonise61(s1); s2=codonise61(s2);
	% disp('s1 and s2 have been codonise61 now')
end

m=size(s1,2);

%select=nargin-2;
select=1;
switch (select)
    case (1)
        disp('optimising everything: t, kappa and omega')
        [para,lnL]=i_optimtko(s1,s2);
        t=para(1);
        kappa=para(2);
        omega=para(3);
        md=modelgy94(omega,kappa);		% build model from optimised values
    case (2)
    	disp('fixed kappa, optimising t and omega.')
        error('under development!')
    case (3)
    	disp('fixed kappa and omega, optimising t.')
        error('under development!')
        kappa=1.6;   % fixed kappa
        omega=0.8;
        md=modelgy94(omega,kappa);
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
stops=TABLE(icode,:)=='*';
TABLE=TABLE(icode,:);
TABLE(stops)=[];
M=zeros(61);
for i=1:61
for j=i:61
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

%%
% Calculate pS and pN when omega = 1
%%
md0=modelgy94(1,kappa);
Q0=composeQ(md0.R,diag(md0.freq))./61;
pS0=sum(sum(Q0.*M));
pN0=1-pS0;

%%
% Calculates dS and dN
%%
dS=t*pS/(pS0*3);
dN=t*pN/(pN0*3);


%%
% Outputs
%%
if (nargout>2),
dN_dS=dN./dS;
value.lnL=lnL;
value.kappa=kappa;
value.omega=omega;

Ss=pS0*3*m;
Ns=3*m-Ss;
value.S=Ss;
value.N=Ns;

	if (noise),
		fprintf('S=%.1f, N=%.1f\n',Ss,Ns);
	end

end

disp(' ')
disp('NOTE: Here we used an equal codon frequence!')
disp('(i.e., in codeml setting: CodonFreq = 0  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table)')
if (kappa>990)
      disp('The estimate of kappa is infinity. This can happen when you have extreme data,')
      disp('with very divergent sequences, or when there is no transversion.');
end
disp(' ')


function [para,lnL]=i_optimtko(s1,s2)
	global noise;
	et=0.5; ek=1.5; eo=0.8;    % initial values for t, kappa and omega
	options = optimset('fminsearch');

	if (noise)
        options=optimset(options,'display','iter');
    else
        options=optimset(options,'display','off');
	end

	[para,f_opt]=fminsearch(@i_likelifuntko,[et,ek,eo],options,s1,s2);
	lnL=-f_opt;
	if (noise),
		fprintf('lnL = %.5f\n',lnL);
		fprintf('t=%.5f, kappa=%.5f, omega=%.5f\n',para(1),para(2),para(3));
	end

function [lnL] = i_likelifuntko(x,s1,s2)
	lnL=inf;
	if (any(x<eps)||any(x>999)), return; end
	t=x(1); kappa=x(2); omega=x(3);
	if (t<eps||t>5), return; end
	if (kappa<eps||kappa>999), return; end
	if (omega<eps||omega>10), return; end
	md=modelgy94(omega,kappa);
	[lnL] = -1*likelidist(t,md,s1,s2);


%
% Fixed omega
%
function [para,lnL]=i_optimtk(s1,s2,omega)
	global noise;
	et=0.5; ek=1.5;    % initial values for t and kappa
	options = optimset('fminsearch');
	if (noise), options=optimset(options,'display','iter'); end
	[para,f_opt]=fminsearch(@i_likelifuntk,[et,ek],options,s1,s2,omega);
	lnL=-f_opt;
	if (noise),
		disp(sprintf('lnL = %.5f',lnL))
		disp(sprintf('t=%.5f, kappa=%.5f, omega(fixed)=%.5f',para(1),para(2),omega))
	end

function [lnL] = i_likelifuntk(x,s1,s2,omega)
	lnL=inf;
	if (any(x<eps)||any(x>999)), return; end
	t=x(1); kappa=x(2);
	if (t<eps||t>5), return; end
	if (kappa<eps||kappa>999), return; end
	md=modelgy94(omega,kappa);
	[lnL] = -1*likelidist(t,md,s1,s2);

%%
% Fixed omega and kappa
%%
function [para,lnL]=i_optimt(s1,s2,kappa,omega)
	global noise;
	et=0.5;  % initial values for t
	options = optimset('fminsearch');
	if (noise), options=optimset(options,'display','iter'); end
	md=modelgy94(omega,kappa);
	[para,f_opt]=fminsearch(@i_likelifunt,[et],options,s1,s2,md);
	lnL=-f_opt;
	if (noise),
		disp(sprintf('lnL = %.5f',lnL))
		disp(sprintf('t=%.5f, kappa(fixed)=%.5f, omega(fixed)=%.5f',para(1),kappa,omega))
	end

function [lnL] = i_likelifunt(x,s1,s2,md)
	lnL=inf;
	t=x(1);
	if (t<eps||t>5), return; end
	[lnL] = -1*likelidist(t,md,s1,s2);
