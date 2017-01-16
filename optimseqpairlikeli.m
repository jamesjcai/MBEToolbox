function [d,lnL] = optimseqpairlikeli(model,s1,s2,ax,bx)
%OPTIMSEQPAIRLIKELI - Optimises distance of given sequence pair under a model
%
% Syntax: [d,lnL] = optimlikelidist(model,s1,s2,ax,bx)
%
% Inputs:
%    model   - Substitution model structure
%    sl      - Sequence 1
%    s2      - Sequence 2
%    ax      - Lower boundary (default = 0)
%    bx      - Upper boundary (default = 2)
%
% Outputs:
%    model.R   - Rate matrix
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


if (nargin<4), ax=0.00001; bx=5; end
if (ax<=eps), ax=0.00001; end

g = inline('-1*seqpairlikeli(t,model,s1,s2)','t','model','s1','s2');
options = optimset('fminbnd');
[toptim,lnLmax]= fminbnd(g,ax,bx,options,model,s1,s2);
d=toptim;
lnL=-lnLmax;