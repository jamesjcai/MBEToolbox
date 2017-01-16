function [t]=translateseq(s,icode)
%TRANSLATESEQ - Translate coding DNA sequence into protein.
%
% Syntax: [t]=translateseq(s,icode)
%
% Inputs:
%    s       - protein-coding nucleotide sequences
%    icode   - Genetic code. 1 - Standard; 2 - Vertebrate Mitochondrial;
%              3 - Yeast Mitochondrial; 4 - Mold, Protozoan, Coelenterate
%              Mitochondrial, and Mycoplasma/Spiroplasma; 5 - Invertebrate
%              Mitochondrial; 6 - Ciliate, Dasycladacean, and Hexamita
%              Nuclear
% Outputs:
%    t     - Protein sequences
%
% See also: TRANSLATEALN

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (nargin<2), icode=1; end
if (isstruct(s)), icode=s.geneticcode; s=s.seq; end

[TABLE,CODON] = codontable;
G=TABLE(icode,:);
t=G(codonise64(s));
t=encodeseq(t,3);

%5'3' Frame 1
%atgaaaacattgaaccttttggcacctgtttcgcaattgtttaagccccttcaaaggcaa
% M  K  T  L  N  L  L  A  P  V  S  Q  L  F  K  P  L  Q  R  Q
%agttcttcacgcgtaaacgcccttttgataaatgttcctagagttcagctttctaaaggc
% S  S  S  R  V  N  A  L  L  I  N  V  P  R  V  Q  L  S  K  G
%aagaacaagtaccttttgatgatgattcacatgcacggtttaactaggtttggcaggacc
% K  N  K  Y  L  L  M  M  I  H  M  H  G  L  T  R  F  G  R  T
%attgtaagaggatccgcctccaaggaccacgaggaaatctttgaggaaatccagaaggag
% I  V  R  G  S  A  S  K  D  H  E  E  I  F  E  E  I  Q  K  E
%atggacaaaattggcatttgcgccaaatgcctgggtggagggttcatcagcaacaaggag
% M  D  K  I  G  I  C  A  K  C  L  G  G  G  F  I  S  N  K  E
%gacaagaaagtaatgaaaatttacgggtgctgcaaaacttttggggaggcgccgcacggc
% D  K  K  V  M  K  I  Y  G  C  C  K  T  F  G  E  A  P  H  G
%aggacgaaggacatactgctgtcttggaccaaattccagcattacaacatcatcttgcc
% R  T  K  D  I  L  L  S  W  T  K  F  Q  H  Y  N  I  I  L