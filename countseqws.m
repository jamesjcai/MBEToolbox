function [ws,sw,ss,ww]=countseqws(s1,s2)
%COUNTSEQWS -
%
% [ws,sw,ss,ww]=countseqws(s1,s2)
%
% s1 is ancestral sequence.
%
% AT->GC changes are referred to as "weak-to-strong" (W->S) as they result
% in a replacement of a weak A:T bond with a "strong" G:C bond.
%
%REF: http://dx.doi.org/10.1371/journal.pbio.1000026

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

TWS = [0, 0, 0, 0;
       1, 0, 0, 1;
       1, 0, 0, 1;
       0, 0, 0, 0];

TSW = [0, 1, 1, 0;
       0, 0, 0, 0;
       0, 0, 0, 0;
       0, 1, 1, 0];
   
TSS = [0, 0, 0, 0;
       0, 0, 1, 0;
       0, 1, 0, 0;
       0, 0, 0, 0];

TWW = [0, 0, 0, 1;
       0, 0, 0, 0;
       0, 0, 0, 0;
       1, 0, 0, 0];
   
[X]=countchange(s1,s2,4);

ws = sum(sum(TWS.*X));
sw = sum(sum(TSW.*X));
ss = sum(sum(TSS.*X));
ww = sum(sum(TWW.*X));

