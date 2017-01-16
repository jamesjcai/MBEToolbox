%% MBEToolbox DEMO - Genetic distances estimation
% Welcome to MBEToolbox.  This is a demonstration of
% MBEToolbox functions for extimating genetic distances.
%
% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


cd(fileparts(which('mbe_demo3')))
% Current working directory will be changed.

%%
% Let's load an example file

filename=fullfile('seq_examples','ng-new.phy');
aln = readphylip_i(filename,2,1)

%%
% Now, calculating Jukes & Cantor 69 distance. Matrix containing
% pairwise distances can be display.

D=dn_jc(aln);
printmatrix(D,aln)

%%
% Kimura 80 Distance

D=dn_k2p(aln);
printmatrix(D,aln)

%%
% Log-det (paralinear) distance

[D]=dn_logdet(aln);
printmatrix(D,aln)

%%
% Log-det (paralinear) distance

alpha=0.5;
[D]=dn_jin_nei90(aln,alpha);
printmatrix(D,aln)

%%
% Compute syn. non-syn. substitutions rates using Nei-Gojobori method

[dS,dN]=dc_ng86(aln);
printmatrix(dS,aln)
printmatrix(dN,aln)

%%
% Compute syn. non-syn. substitutions rates using Li93 method

[dS,dN]=dc_li93(aln);
printmatrix(dS,aln)
printmatrix(dN,aln)


%%
% Matrix containing pairwise distances can be exported into TEXT,
% EXCEL or WORD file.
if (ispc)
   exportdismatrix(D,aln)
else
   disp('This function is only allowed for Windows system.')
end


%%
% This code creates GTR model from scratch, and then uses the model
% to estimate GTR distance and gives likelihood.

s1=[1 2 3 4];   % ACGT
s2=[2 2 3 4];   % CCGT
rmatrix=[1.0 1.33333 1.0 1.0 1.333333 1.0];
freq=[.1 .2 .3 .4];
[model] = modelgtr(rmatrix,freq)
[d,lnL] = optimlikelidist(model,s1,s2,0,5)

%%
% Finally, we show how to use MBEToolbox to do sliding window analysis
% of Ka and Ks.

filename=fullfile('seq_examples','HCV_aligned.fas');
aln = readfasta(filename,2,1);
plotslidingwinkaks(aln);

%%
% Thank you for viewing this introduction to MBEToolbox distances
% estimation funcation.