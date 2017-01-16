function [nk]=basescovered(k,c,G)
%BASESCOVERED - number of bases covered exactly k times
%
%[nk]=basescovered(k,c,G)
%
% nk  - the number of bases covered by exactly k times
% c   - coverage
% G   - genome size
%
% When sequencing a single genome, the Lander–Waterman model based on the
% assumptions of independent and random reads implies that the coverage of
% each base is distributed according to a Poisson distribution with
% parameter c (the coverage). Defining n_k to be the number of bases
% covered exactly k times and G to be the genome size, we have
%
% REF: Bioinformatics for Whole-Genome Shotgun Sequencing of Microbial 
% Communities Kevin Chen*, Lior Pachter*

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if nargin<3
    G=1;
end

nk=G.*(c.^k).*exp(-c)./factorial(k);
nk2=poisspdf(k,c)


% Poisson Calculations for Shotgun Sequencing
%
% A simple formula can be used to estimate what percent of a clone will be
% sequenced for a certain level of random sequencing (Lander & Waterman
% 1988).  The table below shows several examples.  For example, if you
% obtain 100 kb of random sequence from a 100 kb clone (1 x coverage) then
% you would expect to have 63% of  that clone sequenced.  Some regions of
% the clone will have been sequenced several times and other regions will
% not have been sequenced at all. 
%
% Fold coverage 	Percent of clone sequenced
% 0.25 x	22%
% 0.50 x	39%
% 0.75 x	53%
% 1 x	63%
% 2 x	88%
% 3 x	95%
% 4 x	98%
% 5 x	99.4%
% 6 x	99.75%
% 7 x	99.91%
% 8 x	99.97%
% 9 x	99.99%
% 10 x	99.995%
%
% These figures should be used as a guide only and there are many reasons 
% why actual results may deviate from them, perhaps the most significant
% being non-random shotgun cloning.
%
% Lander ES, Waterman MS (1988) Genomic mapping by fingerprinting random
% clones: a mathematical analysis. Genomics, 2, 231-239.

% sum(basescovered([1:100],1))
% 1-poisscdf(0,1)















