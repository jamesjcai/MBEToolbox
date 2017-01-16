function geneplot(params)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (nargin<1)
params=[30417 30944;30974 30995;31039 31377;31450 31605;31900 32453;32514 32642;33329 33334;33381 33404;33633 33674];
end
%params=[2 100; 200 250];
nb = size(params,1);  %  number of boxplots
hw = 0.15;   %  halfwidth of boxes

hold on
for ii = 1:nb
   temp1 = params(ii,1);
   temp2 = params(ii,2);
   xx = [temp1 temp1 temp2 temp2 temp1];
%   yy = [ii-hw ii+hw ii+hw ii-hw ii-hw];
   yy = [1-hw 1+hw 1+hw 1-hw 1-hw];
   %plot(xx,yy,'-')
   %fill(xx,yy,'r');
   patch(xx,yy,'r');
end

plot([min(min(params)) max(max(params))],[1 1],'--')

%  make some extra space
axlim = axis;
axlim(1) = axlim(1)-1;
axlim(2) = axlim(2)+1;

%axlim(1) = 0;
%axlim(2) = max(max(params));

axlim(3) = axlim(3)-2;
axlim(4) = axlim(4)+2;

axis(axlim)
axis off
