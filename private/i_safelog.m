function [y] = i_safelog(x)
x(find(x<=0))=realmin;
y=log(x);