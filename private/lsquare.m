function [b,m,E]=lsquare(x,y)
%inputs: x,y column vectors of same lengths
%outputs intercept b, slope m, and error E in least square fit of data to y=m*x+b
% 
% Description: 
% 
% 1) The function file lsquare.m receives two column vectors of data as input and 
% generates  intercept b, slope m, and error term E of a least square fit of the 
% data to a straight line y=b+m*x as output. Function call: 
% 
% [b,m,E]=lsquare(x,y) 
% 
% In addition data and fitted model line are plotted. Save the file and call it in 
% the workspace or from a script file. 
% 
% 2) The script file run_lsquare.m is the driver file for  lsquare.m. In its present 
% form it solves Problem 5.4(a). 
% http://www.math.colostate.edu/~gerhard/classes/331/lsquare.html

%Computation of b,m
P=length(x);
X=[ones(P,1) x];
vec=X\y;
b=vec(1);
m=vec(2);
%Computation of error
dy=y-m*x-b;
E=dy'*dy;
%Plot results
%plot(x,y,'ko',x,m*x+b,'k-v')
%legend('data','model')
