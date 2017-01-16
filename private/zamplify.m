function [Z1] = zamplify(Z)
[n,m]=size(Z);
Xaxis = 1:m;
Yaxis = Z;
[b,m,E]=lsquare(Xaxis',Yaxis');	% m is slope
Z1=Z-m*Xaxis;