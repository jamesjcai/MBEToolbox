function [s,aln] = cdspairgen(len)
%CDSPAIRGEN - generates a random coding sequence and then mutate it

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (nargin<1)
      len=500;
end


freq_aspergillus_fumigatus=[14.2 23.5 34.11 14.42 11.86 20.98 12.13 13.6 7.37 15.63 5.91 9.52 4.59 28.21 21.08 17.32 14.47 12.9 25.43 11.41 11.29 17.55 12.91 16.36 9.39 15.49 10.83 10.88 6.54 23.43 25.17 15.4 23.09 29.96 36.85 26 14.98 28.87 18.18 22.19 14.32 25.61 10.29 18.7 5.5 25.51 16.84 14.6 18.44 11.44 9.98 18.23 14.9 13.47 7.91 14.71 4.76 4.09 26.13 15.58 13.09];
freq=freq_aspergillus_fumigatus./(sum(freq_aspergillus_fumigatus));

[x]=randseq2(1,len,2,freq);
x=codonise61(x);

omega=0.7; kappa=4; t=2.4;
md = modelgy94(omega,kappa,freq);
[x2]=mutateseq(x,md,t);

x=decodonise61(x);
x2=decodonise61(x2);
s=[x;x2];
if (nargout>1)
	aln.seqtype=2;
	aln.geneticcode=1;
	aln.seqnames={'seq1','seq2'};
	aln.seq=s;
end
