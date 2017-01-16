function [k] = compcomp(s)
%COMPCOMP - Calculates the compositional complexity of sequence
%Compisitional complexity in a sequence of length $L$ is given by (Wootton
%and Federhen 1996)
%$$K = 1/L\log _N(L!/\prod\limits_{{\rm{all~}}i} {n_i!)}$$
%where $N$ is 4 for nucleic acid sequences and 20 for protein sequences, and
%$n_i$ are the numbers of each residue in the sequence. $K$ will vary from
%0 for very low complexity to 1 for high complexity. For the sequence CTGA,
%$L!=4\times3\times2\times1=24$, $n_A=n_C=n_G=n_T=1$, $\prod n_i=1$. All
%$K=1/4\log _4(24/1)=0.573$
%
% for nt seequence:
%
% k = (1/L)log4(L!/(n1!*n2!*n3!*n4!))
%
% Syntax: [k] = compcomp(s)
%
% Inputs:
%    s   - Sequence (nucleotide)
%
% Outputs:
%    k   - Compositional complexity
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


[n,m]=size(s);
k=zeros(n,1);
for (i=1:n)
    ss=s(i,:);
    L=length(ss);
    [N] = ntcomposition(ss);
    % k(i,1) = (1/L)*log4((factorial(L)/i_prodfactorial(N)));
    k(i,1) = (1/L)*log4(multinomial(L,N));
    % by this manipulation we avoid overflows...
end


function y = log4(x)
% Compute y = log2(x)/log2(10) with an averaging process so that roundoff
% errors cancel and log10(10^k) == k for all integers k = -307:308.
% Use y = log2(x)/c1 + log2(x)/c2 where c1 = 2*log2(10)+2.5*eps,
% c2 = 2*log2(10)-1.5*eps are successive floating point numbers on
% either side of 2*log2(10).
y=log(x)/log(4);

function y = log20(x)
y=log(x)/log(20);


function [y] = i_prodfactorial(x)
    len=length(x);
	y=1;
	for (k=1:len),
	      y=y*factorial(x(k));
	end


function [y] = i_product(x)
len=length(x);

if (len==1)
	y=1;
	for (k=1:x),
	      y=y*k;
	end
else
	%y=ones(1,len);
	y=1;
	for (k=1:len),
	      % y(k)=i_product(x(k));
	      y=y*i_product(x(k));
	end

end




function c = multinomial(n,varargin)
% MULTINOMIAL Multinomial coefficients
%
%   MULTINOMIAL(N, K1, K2, ..., Km) where N and Ki are floating point
%   arrays of non-negative integers satisfying N = K1 + K2 + ... + Km,
%   returns the multinomial  coefficient   N!/( K1!* K2! ... *Km!).
%
%   MULTINOMIAL(N, [K1 K2 ... Km]) when Ki's are all scalar, is the
%   same as MULTINOMIAL(N, K1, K2, ..., Km) and runs faster.
%
%   Non-integer input arguments are pre-rounded by FLOOR function.
%
% EXAMPLES:
%    multinomial(8, 2, 6) returns  28
%    binomial(8, 2) returns  28
%
%    multinomial(8, 2, 3, 3)  returns  560
%    multinomial(8, [2, 3, 3])  returns  560
%
%    multinomial([8 10], 2, [6 8]) returns  [28  45]

% Mukhtar Ullah
% November 1, 2004
% mukhtar.ullah@informatik.uni-rostock.de

nIn = nargin;
error(nargchk(2, nIn, nIn))

%if ~isreal(n) || ~isfloat(n) || any(n(:)<0)
%    error('Inputs must be floating point arrays of non-negative reals')
%end

arg2 = varargin;
dim = 2;

if nIn < 3
    k = arg2{1}(:).';
%    if isscalar(k)
%        error('In case of two arguments, the 2nd cannot be scalar')
%    end
else
    [arg2{:},sizk] = sclrexpnd(arg2{:});
    if sizk == 1
        k = [arg2{:}];
    else
        if ~isscalar(n) && ~isequal(sizk,size(n))
            error('Non-scalar arguments must have the same size')
        end
        dim = numel(sizk) + 1;
        k = cat(dim,arg2{:});
    end
end

%if ~isreal(k) || ~isfloat(k) || any(k(:)<0)
%    error('Inputs must be floating point arrays of non-negative reals')
%end

n = floor(n);
k = floor(k);

if any(sum(k,dim)~=n)
    error('Inputs must satisfy N = K1 + K2 ... + Km ')
end

c = floor(exp(gammaln(n+1) - sum(gammaln(k+1),dim)) + .5);



% function y = multcoef(n,x)
% % MULTCOEF Multinomial coefficient.
% %  Y = MULTCOEF(N,X) returns the multinomial coefficient with parameter N at the values in X.
% %  Let {X1, X2,  . . .  , Xk}, k > 1, be a set of random variables, each of which can take
% %  the values 0, 1,  . . .  , n; such that for every set of k nonnegative integers {n1,
% %  . . .  , nk} whose sum is n, the multinomial coefficient is,
% %
% %                                                        (n1 + n2 + ... nk)!
% %        C(n; n1, n2, ..., nk) = (n1, n2, ... nk)! =  -------------------------  .
% %                                                     n1! × n2! ×  . . .  × nk!
% %
% %  It is possible to work with large factorials.
% %
% %  Syntax: y = multcoef(n,x)
% %
% %  Inputs:
% %       n - number of trials.
% %       x - vector of the interested values.
% %  Outputs:
% %       y - multinomial coefficient.
% %
% %  Example. Assume that a die is thrown 60 times (n=60) and a record is kept of the number of
% %  times a 1, 2, 3, 4, 5, or 6 is observed. Each trial (e.g., throw of a die) has a (1 or 2 or
% %  3 or . . . or 6) mutually exclusive outcome. In the die tossing data, k = 6, we are
% %  interested to get the multinomial coefficient.
% %
% %  Calling on Matlab the function:
% %             multcoef(n,x)
% %
% %  where n=60 and x=[13,10,8,10,12,7];
% %
% %  Answer is:
% %
% %  ans
% %      = 1.0425e+042
% %
% %  Created by A. Trujillo-Ortiz, R. Hernandez-Walls and A. Castro-Perez
% %             Facultad de Ciencias Marinas
% %             Universidad Autonoma de Baja California
% %             Apdo. Postal 453
% %             Ensenada, Baja California
% %             Mexico.
% %             atrujo@uabc.mx
% %  Copyright (C) January 23, 2005.
% %
% %  To cite this file, this would be an appropriate format:
% %  Trujillo-Ortiz, A., R. Hernandez-Walls and A. Castro-Perez. (2005). multcoef:
% %    Multinomial coefficient. A MATLAB file. [WWW document]. URL http://
% %    www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6786
% %
% %  References:
% %
% %  Abramowitz, M. and Stegun, I. A. (1964), Handbook of Mathematical
% %           Functions, Government Printing Office, 823.24.1.2. Available on
% %           Internet at the URL address http://hcohl.shell42.com/as/frameindex.htm
% %

% if nargin < 2,
%    error('You need to input two arguments.');
%    return,
% end;

% if (length(n)~=1) | (fix(n) ~= n) | (n < 0),
%    error('n must be a positive integer.');
%    return,
% end;

% if sum(x) ~= n
%     error('Inputs must satisfy n = x1 + x2 ... + xi.');
%     return,
% end;

% factor1 = sum(log(1:n));

% c = length(x);

% f = [];
% for i = 1:c,
%     f = [f log(1:x(i))];
% end;

% factor2 = sum(f);

% y = exp(factor1-factor2);
% y = fix(y);

% return,


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <code>multinomial</code> gives right answer in the following particular case:

% multinomial(4,[1 1 1 1])      = 24 (correct)
% multcoef(4,[1 1 1 1])         = 23 (wrong)

function [varargout] = sclrexpnd(varargin)
% SCLREXPND expands scalars to the size of non-scalars.
%    [X1,X2,...,Xn] = SCLREXPND(X1,X2,...,Xn) expands the scalar
%    arguments, if any, to the (common) size of the non-scalar arguments,
%    if any.
%
%    [X1,X2,...,Xn,SIZ] = SCLREXPND(X1,X2,...,Xn) also returns the
%    resulting common size in SIZ.
%
% Example:
% >> c1 = 1; c2 = rand(2,3); c3 = 5; c4 = rand(2,3);
% >> [c1,c2,c3,c4,sz] = sclrexpnd(c1,c2,c3,c4)
%
% c1 =
%      1     1     1
%      1     1     1
%
% c2 =
%     0.7036    0.1146    0.3654
%     0.4850    0.6649    0.1400
%
% c3 =
%      5     5     5
%      5     5     5
%
% c4 =
%     0.5668    0.6739    0.9616
%     0.8230    0.9994    0.0589
%
% sz =
%      2     3

% Mukhtar Ullah
% November 2, 2004
% mukhtar.ullah@informatik.uni-rostock.de

C = varargin;
issC = cellfun('prodofsize',C) == 1;
if issC
    sz = [1 1];
else
    nsC = C(~issC);
    if ~isscalar(nsC)
        for i = 1:numel(nsC), nsC{i}(:) = 0;  end
        if ~isequal(nsC{:})
            error('non-scalar arguments must have the same size')
        end
    end
    s = find(issC);
    sz = size(nsC{1});
    for i = 1:numel(s)
        C{s(i)} = C{s(i)}(ones(sz));
    end
end
varargout = [C {sz}];
