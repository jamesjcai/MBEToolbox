function S=i_encode_a(Seq)
%I_ENCODE_A - Convert an amino acid sequence from a letter to an integer representation
%Internal function encodes the nucleotide sequences by digits
%
% Syntax: S=i_encode_n(Seq)
%
% Inputs:
%    Seq   - Letter representation of sequence
%
% Outputs:
%    S     - Integer representation of sequence
%
%
% See also: I_ENCODE_N 

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

[n,m]=size(Seq);
S = zeros(n,m);  i = 1; j = 1;

code = zeros(256,1); 
for o = 1:256
    code(o) = 22;	% 22 = '-' gaps and any other characters
end

[NT,AA] = seqcode;
for o = 1:21		% 21 = '*' stop codon
   code(abs(AA(o))) = o;  
end

for i=1:n
	for j=1:m
		S(i,j) = code(abs(Seq(i,j)));
	end;
end;
% S=uint8(S);