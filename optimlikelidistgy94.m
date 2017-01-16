function [t,omega,kappa,lnL]=optimlikelidistgy94(s1,s2)
%OPTIMLIKELIDISTGY94 - Optimises distance (t) under a GY94 codon model
%
% Syntax: [t,omega,kappa,lnL]=optimlikelidistgy94(s1,s2)
%
% Inputs:
%    sl     - Sequence 1
%    s2     - Sequence 2
%
% Outputs:
%    t       - Optimised distance
%    lnL     - Maximum log-likelihood
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


%[lnL,X,Y,Z] = optimlikelidistk2p(s1,s2,ax,bx,k1,k2)
%optimlikelidist

%if (nargin<3), flag=0; end


%omega=1;
estKappa = estimatekappa([s1;s2],'hky');

[s1] = codonise61(s1);
[s2] = codonise61(s2);


options = optimset('fminsearch');
[xy,f_opt]=fminsearch(@i_likelidist_gy94,[0.5,0.2,estKappa],options,s1,s2);
t=xy(1);
omega=xy(2);
kappa=xy(3);
lnL=-f_opt;


function [lnL] = i_likelidist_gy94(x,s1,s2)
     t=x(1);
     kappa=x(3);
     omega=x(2);
     if (t<=0|kappa<=0)
	lnL=inf; return;
     end
     if (omega<=0|omega>=10)
	lnL=inf; return;
     end
     [model] = modelgy94(omega,kappa);
     [lnL] = -1*likelidist(t,model,s1,s2);


% [xy,f]=fminsearch('fBMEN289_9',[2.5,1])
% [xy,f_opt]=fminsearch(@fBMEN289_9,[2.5,1])
% [xy,f]=fminsearch(inline('-(2*x(1)*x(2)+2*x(1)-x(1)^2-2*x(2)^2)'),[5,1])
%
%function f=fBMEN289_9(x)
%f=-(2*x(1)*x(2)+2*x(1)-x(1)^2-2*x(2)^2);
%format long


function [kappa] = i_estimatekappa(S)
	[p,q]=countseqpq(S);
	[n,m] = size(S);

	p=p./m; q=q./m;
	w1 = 1-2*p-q;
	w2 = 1-2*q;
	r=(-log(w1)./-log(w2))-0.5;
	kappa=r*2;


%NOTE:
% > I find that results from fminsearch are: 1)highly dependent on
% > initial position in parameter space and 2)not always terminated at
% > the global minimum (known from a grid search). I have played around
% > with the tolerances to no avail.
% >
% > Have others had similar experiences with fminsearch? Any suggestions?
%
% This is what makes global minimization hard to solve.
%
% Imagine you are walking a terrain yourself and you follow
% the terrain down to a deep valley. How do you know you
% didn't miss a deeper valley on the other side of that
% ridge 20 km west? Now make the problem exponentially worse
% by going to n search dimensions instead of 2.
%
% One strategy to improve your *chances* of finding the
% global optimum is to introduce some carefully designed
% element of randomness into the search, so that with
% some probability you will in fact go past that
% ridge. Simulated annealing and genetic algorithms are
% two such strategies.