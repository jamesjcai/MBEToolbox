function out = slidingmean(in, N, step)

if (isempty(in)) | (N<=0)                                              % If the input array is empty or N is non-positive,
 disp(sprintf('SlidingAvg: (Error) empty input data or N null.'));     % an error is reported to the standard output and the
 return;                                                               % execution of the routine is stopped.
end % if

if (N==1)                                                              % If the number of neighbouring points over which the sliding 
 out = in;                                                             % average will be performed is '1', then no average actually occur and
 return;                                                               % OUTPUT_ARRAY will be the copy of INPUT_ARRAY and the execution of the routine
end % if                                                               % is stopped.

nx   = length(in);             % The length of the input data structure is acquired to later evaluate the 'mean' over the appropriate boundaries.

if (N>=(2*(nx-1)))                                                     % If the number of neighbouring points over which the sliding 
 out = mean(in)*ones(size(in));                                        % average will be performed is large enough, then the average actually covers all the points
 return;                                                               % of INPUT_ARRAY, for each index of OUTPUT_ARRAY and some CPU time can be gained by such an approach.
end % if                                                               % The execution of the routine is stopped.



out = zeros(size(in));         % In all the other situations, the initialization of the output data structure is performed.

if rem(N,2)~=1                 % When N is even, then we proceed in taking the half of it:
 m = N/2;                      % m = N     /  2.
else                           % Otherwise (N >= 3, N odd), N-1 is even ( N-1 >= 2) and we proceed taking the half of it:
 m = (N-1)/2;                  % m = (N-1) /  2.
end % if

for i=1:nx,                                                 % For each element (i-th) contained in the input numerical array, a check must be performed:
 if ((i-m) < 1) & ((i+m) <= nx)                             % If not enough points are available on the left of the i-th element..
  out(i) = mean(in(1:i+m));                                 % then we proceed to evaluate the mean from the first element to the (i + m)-th.
 elseif ((i-m) >= 1) & ((i+m) <= nx)                        % If enough points are available on the left and on the right of the i-th element..
  out(i) = mean(in(i-m:i+m));                               % then we proceed to evaluate the mean on 2*m elements centered on the i-th position.
elseif ((i-m) >= 1) & ((i+m) > nx)                          % If not enough points are available on the rigth of the i-th element..
  out(i) = mean(in(i-m:nx));                                % then we proceed to evaluate the mean from the element (i - m)-th to the last one.
 elseif ((i-m) < 1) & ((i+m) > nx)                          % If not enough points are available on the left and on the rigth of the i-th element..
  out(i) = mean(in(1:nx));                                  % then we proceed to evaluate the mean from the first element to the last.
 end % if
end % for i
