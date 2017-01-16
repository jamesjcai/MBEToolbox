function [pt] = bary2pt(coord)
%BARY2PT - Compute the point, pt from barycentric coordinates

%Barycentric coordinates are triples of numbers (t1,t2,t3) corresponding to 
%masses placed at the vertices of a reference triangle (A1,A2,A3). These masses 
%then determine a point P, which is the geometric centroid of the three masses, 
%and is identified with coordinates (t1,t2,t3).
%
%REF: http://mathworld.wolfram.com/BarycentricCoordinates.html

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

coord = coord./sum(coord);
a=[0 0.5 1]; b=[0 1 0];
ab=[a;b]';
pt=zeros(1,2);
for (k=1:3),
      pt=pt+ab(k,:)*coord(k);
end
