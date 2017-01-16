function S=i_encode_n(Seq)
%I_ENCODE_N - Convert a nucleotide sequence from a letter to an integer representation
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
% See also: I_ENCODE_A

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

method=1;
[NT,AA] = seqcode;
switch (method)
    case (1)
	[n,m]=size(Seq);
	S = ones(n,m).*5;
	Seq(Seq=='U')='T';    %replace U with T
	for k=1:5
		S(Seq==NT(k))=k;
    end
    
    case (2)
	[n,m]=size(Seq);
	S = zeros(n,m);  i = 1; j = 1;

	code = zeros(256,1); 
	for o = 1:256
	    code(o) = nan;
	end

	% NT = 'ACGT-';
	for o = 1:5
	   code(abs(NT(o))) = o;  
	end

	for i=1:n
		for j=1:m
			if Seq(i,j) == 'U' 
				S(i,j) = code(abs('T'));
			else
				S(i,j) = code(abs(Seq(i,j)));
			end;
		end;
	end;
end
