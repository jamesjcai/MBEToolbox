%% MBEToolbox DEMO - Phylogenetic inferences
% Welcome to MBEToolbox.  This is a demonstration of
% MBEToolbox functions for phylogenetic inferences.
%
% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



cd(fileparts(which('mbe_demo4')))
% Current working directory will be changed.

%%
% Let's load an example file

filename=fullfile('seq_examples','ng-new.phy');
aln = readphylip_i(filename,2,1)

%%
% Let's performs neighbor joining (NJ) on distance data and plot
% NJ tree. Here we use Log-Det distance.

D=dn_logdet(aln);
printmatrix(D,aln)
plotnjtree(D,1,char(aln.seqnames));

%%
% Build a REV nucleotide substitution model

rmatrix = [1.0 1.33333 1.0 1.0 1.33333 1.0];
freq = [.1 .4 .2 .3];
[model]=modelrev(rmatrix,freq);

%%
% Plot log-likelihood as the function of distance between two sequences

s1='CCAT'; s1=i_encode_n(s1);
s2='CCGT'; s2=i_encode_n(s2);
t=0:0.01:0.8;
lnL=zeros(1,length(t));
for (k=1:length(t)),
      [lnL(k)]=likelidist(t(k),model,s1,s2);
end
plot(t,lnL);
grid;
xlabel('Distance (substitutions/site)');
ylabel('ln(Likelihood)')
hold

%%
% Optimise distance so that we can get maximum likelihood

[optim_t,max_likli] = optimlikelidist(model,s1,s2)
plot(optim_t, max_likli,'r*')

max_likli2 = likelidist(optim_t,model,s1,s2);
max_likli==max_likli2

%%
% MBEToolbox provides interfaces of PHYLIP programs, dnapars, dnaml,
% protpars and proml. Firstly, let's try dnapars...

runphylip(aln,'dnapars')

%%
% ... after translating into protein sequences, let's try proml...
% This may take a while ...

[aln2] = translatealn(aln)
runphylip(aln2,'proml')

%%
% Thank you for viewing this introduction to MBEToolbox phylogenetic
% inference funcations.