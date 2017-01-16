%% MBEToolbox DEMO - New features in version 2.0
% Welcome to MBEToolbox.  This is a demonstration of
% new MBEToolbox functions in vesroin 2.0.
%
% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://www.hku.hk/jamescai/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



cd(fileparts(which('mbe_demo5')))
% Current working directory will be changed.

%%
% Let's load an example protein alignment and tree, and make a model.

if (isunix), slashid='/'; else slashid='\'; end
matfilename=['seq_examples',slashid,'rateofsites_example'];
load(matfilename)

model=modeljtt;

%%
% Let's performs site-specific evolutioanry rate inference by maximum likelihood (ML)
% method.

h = waitbar(0,'Please wait...'); waitbar(30/100,h)
[rateraw,ratenorm]=rateofsites_ml(seq,tree,modeljtt);
close(h)

%%
% We plot estimated rates.

bar(ratenorm)
ylabel('Normalised Relative Evolutionary Rate');
xlabel('Site Position');
axis([0 length(ratenorm)+1 min(ratenorm)*1.1 max(ratenorm)*1.1]);


%%
% Now we demostrate one of recombination detecting method - reticulate.

filename=['seq_examples',slashid,'ng-new.phy'];
aln = readphylip_i(filename,2,1)

reticulate(aln);

%%
% Maximum likelihood estimation of dS and dN. This funciton is slow so we only
% do it for one pair.

h = waitbar(0,'Please wait...'); waitbar(30/100,h)
[dS,dN,dN_dS]=dc_ml(aln,2,9);
close(h)
disp(' ')
disp('===================================')
disp('        dS        dN     dN/dS')
fprintf('%10.4f%10.4f%10.4f\n', dS, dN, dN_dS);
disp('====================================')
disp(' ')

%%
% Thank you for viewing this introduction to MBEToolbox phylogenetic
% inference funcations.