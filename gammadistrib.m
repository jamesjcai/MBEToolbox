function [gm]=gammadistrib(n,a,se)
%GAMMADISTRIB - Model of rate heterogeneity: Discrete Gamma
%  DiscreteGamma
%   Discretization of gamma distribution with equal proportions in each
%   category.
%construct discrete Gamma distribution (mean = 1.0)
%discrete Gamma distribution (Z. Yang. 1994. JME 39:306-314)
% /* discretization of gamma distribution with equal proportions in each
%   category
%
% Syntax: [gm]=gammadistrib(n,a,se)
%
% Inputs:
%    n         - Number of rate categories
%    a(lpha)   - Shape parameter
%    se        - S.E.
%
% Outputs:
%    gm   - Discrete Gamma model
%
% See also: INVDISTRIB

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if(nargin<1), n=4; end

if (nargin<2), a=0.5; end

showSE=0;
if (nargin==3)
  showSE=1;
end


% number of rate categories
gm.numRates=n;
gm.alpha=a;


% probability of each rate
rate=zeros(1,n);
border=zeros(1,n);
for k=1:n
  rate(1,k) = gaminv((2.0*(k-1)+1.0)/(2.0*n), a, 1/a);
  border(1,k) = gaminv(k/n,a,1/a);
end

rate=rate./mean(rate);
gm.rate=rate;

gm.prob=ones(1,n)./n;

if (showSE)
	gm.alphaSE = se;
end


if (showSE)
    fprintf('Model of rate heterogeneity: Discrete Gamma\n');
    fprintf('Number of rate categories: %d\n', gm.numRates);
    fprintf('Gamma distribution parameter alpha: %2.2f\n', gm.alpha);
	fprintf(' (S.E. %2.2f)\n', gm.alphaSE);
    xmin=0; xmax=2;
    
    i_gammaplot2(xmin,xmax,a,1/a);
    hold on;
    for k=1:n-1
        x=border(1,k);
        y=gampdf(x,a,1/a);
        stem(x,y,'r.-');
    end
    %	axis([0 2 0 4]);
    leg=sprintf('alpha=%2.2f',gm.alpha);
    legend(leg)
    title(sprintf('Discrete gamma distribution, n=%d',gm.numRates))
    hold off
end

%getNumParameters = 1;
%getLowerLimit=0.001
%getUpperLimit=100.0
%getDefaultValue=0.5




function i_gammaplot(xmin,xmax,a,b,colr)
%gammaplot(xmin,xmax,a,b,col)
%plots a GAMMA distribution with parameter a and b > 0, in the interval [xmin,xmax]
%col is a text string which indicates graph colour (default col='b')
%Copyright gianlucabaio2002

if nargin==4
   colr='b';
end

x=[xmin:.01:xmax];
%const1=b^a;
%const2=gamma(a);
for i=1:length(x)
%   y(i)=const1/const2*x(i)^(a-1)*exp(-x(i)*b);
%   y(i)=(b^a)/gamma(a)*x(i)^(a-1)*exp(-x(i)*b);
    y(i)=gammad(x(i),a,b);
end
plot(x,y,colr)


function y=gammad(x,alpha,beta)
%GAMMAD	y=GAMMAD(x,alpha,beta) is the density of a Gamma(alpha,beta)
%	random variable at x.

%	GKS 15 May 92
y=x.^(alpha-1).*exp(-x./beta)./gamma(alpha)./(beta.^alpha);



function i_gammaplot2(xmin,xmax,a,b,colr)
%gammaplot(xmin,xmax,a,b,col)
%plots a GAMMA distribution with parameter a and b > 0, in the interval [xmin,xmax]
%col is a text string which indicates graph colour (default col='b')
%Copyright gianlucabaio2002
if nargin==4
   colr='b';
end

x=xmin:.01:xmax;
y=zeros(1,length(x));
for i=1:length(x)
   y(i)=gampdf(x(i),a,b);
end
plot(x,y,colr)
