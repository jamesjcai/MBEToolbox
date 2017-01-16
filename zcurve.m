function [X,Y,Z,Z1,AT,GC]=zcurve(Seq)
%ZCURVE - Return Z curve components
%REFERENCE:
% C.T. Zhang & R. Zhang, (1991) Analysis of distribution of bases in the
%      coding sequences by a diagrammatic technique.
%      Nucl. Acids Res., 6313-6317.
% R. Zhang & C.T. Zhang, (1994) Z Curves, an Intuitive Tool for Visualizing
%      and Analyzing DNA sequences.
%      J. Biomol. Struc. Dynamics 11, 767-782.
%
% Syntax: [X,Y,Z,Z1,AT,GC]=zcurve(Seq)
%
% Inputs:
%    Seq     - Nucleotide sequence, i.e., a vector representation of DNA sequence
%
% Outputs:
%    X    - The x-component of a Z curve, which displays the distribution of
% 	    purine/pyrimidine (R/Y) bases along the sequence
%    Y    - The y-component of a Z curve  which displays the distribution of
% 	    amino/keto (M/K) bases along the sequence
%    Z    - The z-component of a Z curve  which displays the distribution of
% 	    strong-H bond/weak-H bond (S/W) bases along the sequence
%    Z1   - G+C content distribution along the sequence
%    AT   - (X+Y)/2: AT disparity (Xn+Yn)/2 is called the AT-disparity, which represents
% 	    the excess of A over T in the subsequence constituted from the 1st
% 	    base to the nth bases of the sequence being studied.
%    GC   - (X-Y)/2: GC disparity (Xn-Yn)/2 is called the GC-disparity, which represents
%	    the excess of G over C in the subsequence constituted from the 1st
%	    base to the nth bases of the sequence being studied.
%
% See also: PLOTZCURVE

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[n,m] = size(Seq);
X=zeros(1,m); Y=zeros(1,m); Z=zeros(1,m);
Z1=zeros(1,m); AT=zeros(1,m); GC=zeros(1,m);

% NT = 'ACGT-';
for (k=1:m)
	s = Seq(1,1:k);
	X(1,k)= 2*(sum(s==1)+sum(s==3))-k;
	Y(1,k)= 2*(sum(s==1)+sum(s==2))-k;
	Z(1,k)= 2*(sum(s==1)+sum(s==4))-k;
end

Xaxis = 1:m;
Yaxis = Z;
[b,m,E]=lsquare(Xaxis',Yaxis');	% m is slope

Z1=Z-m*Xaxis;
AT=(X+Y)./2;
GC=(X-Y)./2;