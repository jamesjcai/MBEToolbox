function [Y,I] = randmat(X,DIM)
%RANDMAT    Shuffles array elements or cell elements randomly.
%   For vectors, RANDMAT(X) shuffles the array by calling randperm - use
%   randperm instead.
%   For matrices, RANDMAT(X) shuffles the elements randomly.
%   RANDMAT works for 1 or 2 dimensions.
%
%   Y = RANDMAT(X,DIM)
%   has an optional parameter.
%   DIM selects a dimension to shuffle randomly.
%   The result is in Y which has same shape and type as X.
%
%   [Y,I] = RANDMAT(X,DIM) also returns an index matrix I.
%   If X is a vector, then Y = X(I).
%   If X is an m-by-n matrix and DIM=1, then
%       for j = 1:n, Y(:,j) = X(I(:,j),j); end
%
%
%   Example: if X = [1 2 3
%                    4 5 6]
%
%   then RANDMAT(X,1) might be [1 5 3 and RANDMAT(X,2) might be [3 2 1
%                                4 2 6]                          4 6 5]
%
%   See also PERMS, RAND.

%  Author: Peter Bodin <pbodin@kth.se>
%
%  2005-03-21 22:06:44 Peter Bodin <pbodin@kth.se>
%  * Initial revision

% BEGIN LICENSE
%    Copyright (C) 2005 Peter Bodin
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA, 02111-1307 USA
% END LCENSE

msg = nargchk(1, 2, nargin);                            % check input
error(msg);
if ndims(X)>2
    error('RANDMAT works for 1 or 2 dimensions only.')
end
rand('state',sum(100*clock));                           % set rand state
[r,c] = size(X);                                        % size
switch nargin                                           % argument switch
    case 1 || isvector(X)
        I = reshape(randperm(numel(X)),[r c]);          % use randperm if vector
        Y = X(I);                                       % or if all elements should be shuffled
    case 2
        switch DIM                                      % dimension switch
            case 1
                [I,I] = sort(rand([r c]),1);            % random row permutations
                Y = X(I+ repmat([0 (r+1:r:c*r)-1],r,1));
            case 2
                [Y,I] = randmat(X',1);                  % use transposed input for columns
                Y = Y';
            otherwise
                error('RANDMAT works for 1 or 2 dimensions only.')
        end
end

