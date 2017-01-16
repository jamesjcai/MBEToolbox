function [d,F] = uni_dist(s1,s2)
[D]=countBaseChange(s1,s2);

test=0;
if (test)
D=[  261     1     9      2;
      1   267     0     29;
     14     0    82      0;
      1    23     0    206];
end




D=D./sum(sum(D));
for (i=1:4),
for (j=i:4),
if (i~=j)
	D(i,j)=(D(i,j)+D(j,i))./2;
	D(j,i)=D(i,j);
end
end
end




x=trace(D)./4;
y=(sum(sum(triu(D,1)))+sum(sum(tril(D,-1))))./12;
F=zeros(4);
for (i=1:4),
for (j=i:4),
	if (i~=j)
		F(i,j)=y;
		F(j,i)=F(i,j);
	else
		F(i,j)=x;
	end
end
end
%F=F./sum(sum(F));
F
sum(sum(F))




freq=[.25 .25 .25 .25];
PI=diag(freq);
%d=-trace(PI*logm((PI^-1)*D));
d=-trace(PI*logm(inv(PI)*F));




X=D;
model=modeljc;
F=((sum(sum(X))-trace(X))*model.R)./4;
F=eye(4)*trace(X)./4+F;
F
sum(sum(F))
X=-trace(PI*logm(inv(PI)*F))
