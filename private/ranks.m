function[F] = ranks(X,m,ties)
% Purpose  : Rank non-missing elements of a vector into m fractiles
% Inputs   : X    - (n x 1) or (1 x n) vector, possibly including NaNs and Infs
%            m    - number of fractiles; set Inf to produce raw ranks 
%            ties - if 1, ties are given equal raw ranks, set at group minimum
%                   if 0, ties are given consecutive raw ranks, hence fractile
%                         assignments may differ
% Output   : F - (n x 1) or (1 x n) vector of fractile assignments
% Notes    : 1. Fractile 1 collects smallest values of X.
%            2. NaNs are given NaN category assignments and ignored in calculations.
%            3. If r is Xi's raw rank in the non-missing part of X, of length p, 
%               Xi is assigned to fractile f = 1 + floor[r*m/(p+1)].
% Example  : x = rand(5,1); [x ranks(x,5,0)]
%            ans =
%                   0.4204    2.0000
%                   0.8920    4.0000
%                   0.3457    1.0000
%                   0.8084    3.0000
%                   0.9943    5.0000
% Author   : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 3/13/05
% See also : FRACTILE (on Mathworks File Exchange)
s = size(X);
if s(1) > 1 & s(2) > 1
   error('X must be a vector')
end
if ties ~= 0 & ties ~= 1
   error('Parameter Ties must be set to 1 (equal raw ranks) or 0 (consecutive raw ranks)')
end
n = length(X);
Z = X(~isnan(X)); 
p = length(Z);
if ~isinf(m) & length(unique(Z)) < m
   warning('Number of categories exceeeds number of distinct non-missing values')
else
   if ~isinf(m) & p < m
      warning('Number of categories exceeds number of non-missing values')
   end   
end   
Y = sort(X);          
I = 1:length(Y);
R = X;   % raw ranks
for i = 1:n
    Ii = I(X(i) == Y);
    if ~isempty(Ii) 
       j = Ii(1); 
       R(i) = j;
       if ties
          Y(Ii(1)) = NaN;
       end   
    else
       R(i) = NaN;
    end
end   
if isinf(m), F = R;
else         F = 1+floor(R*m/(p+1));
end   
