function [model] = modelq3(rmatrix,freq)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (nargin<2),
    rmatrix=[2 3 1 4 5 6];
    freq=[.3 .2 .4];
    disp('Using example parameters, rmatrix and freq.')
end

model.name='q3 (non-reversible)';
model.R=i_r2R(rmatrix);
if sum(freq(:))-1>eps
    error('freq error')
end
model.freq=freq;


    
    

    

%%
    
function [R] = i_r2R(r)

% r=[a,b,c,d,e,1]
%
% R=
% x a b
% d x c
% e f x 

if length(r)==5
    r=[r,1];
elseif length(r)==6
    r=r./r(6);
else
    error('Wrong r vector')
end

map=[1,2; 1,3; 2,3; 2,1; 3,1; 3,2];

R=zeros(3);
for k=1:6
    id=map(k,:);
    R(id(1),id(2))=r(k);
end
R=3*R./sum(R(:));   % normalized rate


