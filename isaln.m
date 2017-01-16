function [y]=isaln(aln)
%ISALN - Is an alignment structure?
%
% [y]=isaln(aln)

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



y=false;
if ~(isstruct(aln)), return; end
if ~(isfield(aln,'seqtype')), return; end
if ~(isfield(aln,'seqnames')), return; end
if ~(isfield(aln,'seq')), return; end
if ~(isfield(aln,'geneticcode')), return; end

if ~(iscell(aln.seqnames)), return; end
if ~(isnumeric(aln.seqtype)), return; end
if ~(isnumeric(aln.seq)), return; end
if ~(isnumeric(aln.geneticcode)), return; end

if ~(ismember(aln.seqtype, [1 2 3])), return; end
if ~(ismember(aln.geneticcode, 1:13)), return; end
y=true;
