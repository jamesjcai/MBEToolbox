function mutcube(s)

s=codonise61(s)

[n,m]=size(s);
pickn=nchoosek(n,3);         % number of triples
pickv=nchoosek([1:n],3);     % pickers of triples
figure;
for (k=1:pickn)
  x=s(pickv(k,:),:);         % triple seq;
  plot3(x(1,:),x(2,:),x(3,:),'o')
  hold on
end


pickm=nchoosek(m,3);         % number of triples
pickv=nchoosek([1:m],3);     % pickers of triples
figure;
for (k=1:pickm)
  x=s(:,pickv(k,:));         % triple seq;
  plot3(x(:,1),x(:,2),x(:,3),'ro')
  hold on
end




