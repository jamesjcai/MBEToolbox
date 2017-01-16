function [d] = mcmcjc69(theta0)
%MCMCJC69 - Example of distance estimation under jc69 model to explain MCMC algorithm
%
%REF: Mathematics of Evolution and Phylogen, Edited by Olivier Gascuel

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (nargin<1), 
	theta0=.5;  % initial state
end

%omega=.01;
omega=.2;

n=100;
x=10;
u=1/.1;

nstep=5000;
d=zeros(1,nstep);

way=2;
for k=1:nstep
	theta=theta0-omega/2+omega*rand;
	if (theta<0), theta=-1*theta; end
	switch (way)
	    case (1)
            p1=prior(theta,u)/prior(theta0,u);
            p2=jclike(theta,n,x)/jclike(theta0,n,x);
            alpha=min(1, p1*p2);
	    case (2)
            p1=log( prior(theta,u)*jclike(theta,n,x));
            p2=log( prior(theta0,u)*jclike(theta0,n,x));
            alpha=min(1, exp(p1-p2));
	end
	if (rand<alpha), theta0=theta;	end
	d(k)=theta0;
end

subplot(3,1,1),plot(d(1:500))
subplot(3,1,2),plot(d)
subplot(3,1,3),hist(d(500:end),200)


function [p] = prior(theta,u)
	p=u*exp(-u*theta);

function [lik] = jclike(theta,n,x)
	p=.75*(1-exp((-4/3)*theta));
	%lik=nchoosek(n,x)*(p^x)*(1-p)^(n-x);
	lik=(p^x)*(1-p)^(n-x);
    
    