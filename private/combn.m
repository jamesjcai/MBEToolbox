function [M,IND] = combn(V,N)
% COMBN - all combinations of N elements taken from the vector V
% M = COMBN(V,N) returns all combinations of N elements of the elements in
% vector V. M has the size (length(V).^N)-by-N.
%
% [M,I] = COMBN(V,N) also returns the index matrix I so that M = V(I).
%
% V can be an array of numbers, cells or strings.
%
% Example:
%   M = COMBN([0 1],3) returns the 8-by-3 matrix:
%     0     0     0
%     0     0     1
%     0     1     0
%     0     1     1
%     ...
%     1     1     1
%
% All elements in V are regarded as unique, so M = COMBN([2 2],3) returns 
% a 8-by-3 matrix with all elements equal to 2.
%
% NB Matrix sizes increases exponentially at rate (n^N)*N.
%
% See also PERMS

% Jos, 2005

if nargin ~=2,
    error('Two arguments required.') ;
end

if ~isnumeric(N) | numel(N) ~= 1 | N < 1 | fix(N) ~= N,
    error('Second argument should be a positive integer.') ;    
end

if N==1,    
    % special case in which ndgrid is of no use, just copy the vector
    M = V(:) ;
    IND = [1:length(M)]' ;
else
    nV = numel(V) ;
    % don waste space, if only one output is requested
    [M{1:N}] = ndgrid(1:nV) ;
    M = fliplr(reshape(cat(ndims(M{1}),M{:}),[],N)) ;
    if nargout==2,
        IND = M ;
        M = V(IND) ;
    else
        M = V(M) ;
    end
end
