function p = normcdf(x,mu,sigma)

if nargin < 2, mu = 0; end
if nargin < 3, sigma = 1; end

[errorcode x mu sigma] = distchck(3,x,mu,sigma);
if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

p = 0.5 * erfc(-(x-mu)./(sqrt(2)*sigma));
