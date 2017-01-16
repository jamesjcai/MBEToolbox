function [D]=dp_krkc(aln)
%DP_KRKC - Radical/Conservative Replacement Rate Ratio
% Radical/Conservative Replacement Rate Ratio (KR/KC) depends on the
% definition of radical and conservative changes in the classification of amino acids
% REF:  Kousuke Hanada, Shin-Han Shiu, and Wen-Hsiung Li,
%       Mol Biol Evol 2007 24: 2235-2241; doi:10.1093/molbev/msm152
%
% Zhang's method
% http://www.springerlink.com/content/869dxh45fndjekgl/
%
% Modified Zhang's method
% http://www.springerlink.com/content/1d24jrjp3yr41p11/
%
% Syntax: [D]=dp_krkc(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    D     - Distance matrix
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $




function [seqdeg] = i_degseq(seqrel,deglevel)
seqdeg=seqrel;
if (nargin<2), deglevel=1; end
switch (deglevel)
    case (1)
        % Classification (A) by the maximum correlation with the KA/KS
        % ratio
	setaa={encodeseq('ANCGPST',3), encodeseq('ILMV',3),...
        encodeseq('RQHKFWY',3), encodeseq('DE',3) };
    case (2)
        % Classification (B) by polarity and volume
  	  setaa={encodeseq('C',3), encodeseq('AGPST',3), encodeseq('NDQE',3),...
      encodeseq('RHK',3), encodeseq('ILMV',3), encodeseq('FWY',3)};
    case(3)
        % Classification (C) by charge and aromatic
	setaa={encodeseq('DE',3), encodeseq('QAVLICSTNGPM',3),...
        encodeseq('FYW',3),encodeseq('KRH',3)};
    case(4)
        % Classification (D) by charge and polarity
	setaa={encodeseq('STYCNQ',3), encodeseq('DE',3),...
        encodeseq('KRH',3),encodeseq('GAVLIFPMW',3)};

end
for (k=1:length(setaa)),
      x=setaa{k};
      x1=x(1);
      seqdeg(ismember(seqdeg,x))=x1;
end
