function [coord] = pt2bary(pt)
%PT2BARY - Compute barycentric coordinates of the point, pt.
%Compute barycentric coordinates of the point, pt, relative to the three vertex points.

%Barycentric coordinates are triples of numbers (t1,t2,t3) corresponding to 
%masses placed at the vertices of a reference triangle (A1,A2,A3). These masses 
%then determine a point P, which is the geometric centroid of the three masses, 
%and is identified with coordinates (t1,t2,t3).
%REF: http://mathworld.wolfram.com/BarycentricCoordinates.html

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

 
a=[0 0.5 1]; b=[0 1 0];


x = pt(1) ;
y = pt(2) ;
 
A = [ 1, 1, 1 ; a ; b ] ;  % Set up the matrix for the system.
 
coord = A\[1 ; x ; y ] ;   % Solve the linear system.
coord = coord';



%   This function determines whether the point with
%   coordinates pt is inside or outside of the triangle
%   given by the three vertex points with coordinates
%   (a(1),b(1)), (a(2),b(2)), and (a(3),b(3)).
%   The method used is to compute the barycentric coordinates
%   of the point relative to the three vertex points.  The array
%   coord contains the barycentric coordinates.  If all of the 
%   barycentric coordinates are nonnegative, then the point is
%   said to be inside the triangle.
%   If in_out is 1, the point is inside,
%   if in_out is 0, the point is outside,and
%   if in_out is 0.5, the point is on the edge.

% if all( coord > 0 )     % If all barycentric coordinates are
%     in_out = 1 ;        % non-negative, the point is inside.
% elseif any( coord < 0 ) % If any coordinate is negative,
%     in_out = 0 ;        % the point is outside.
% else                    % Otherwise,
%     in_out = 0.5 ;      % it is on the boundary.
% end
