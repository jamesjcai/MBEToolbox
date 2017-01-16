function T = freqtable(x)
%FREQTABLE Frequency table.
%   TABLE = FREQTABLE(X) takes a vector X and returns a matrix, TABLE.
%   The first column of TABLE contains the unique values of X.  The
%   second is the number of instances of each value. 
%
%Note:  Since FREQTABLE employs only built-in functions, it is much faster 
%  than both TABULATE(Statistics Toolbox) and the following one line code:
%
%  [table(:,2), table(:,1)] = hist(x(:), unique(x));

% Mukhtar Ullah
% December 28, 2004
% mukhtar.ullah@informatik.uni-rostock.de

% Rewritten as suggested by urs (us) schwarz (us@neurol.unizh.ch).

%------------------------------------------------------------------------
% To appreciate the cleverness of urs schwarz, I am including below my 
% previous code for interested ones.
%------------------------------------------------------------------------
% x = sort(x(:));
% dx = diff(x);
% idx = find(~dx);
% if isempty(idx)
%     T = [x ones(size(x))];
% else
%     x1 = x([dx; 1] & [1; dx]);
%     k = find(diff(idx) > 1);
%     i2 = idx([1; k+1]);
%     T = [x1 ones(size(x1)); x(i2) idx([k; end]) - i2 + 2];
%     [T(:,1), i1] = sort(T(:,1));
%     T(:,2) = T(i1,2);
% end
%------------------------------------------------------------------------

x = sort(x(:));
eg = x([find(diff(x));end]);
T = [eg histc(x,eg)];