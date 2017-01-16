function writepaml(aln,filename)
%WRITEPAML - Write alignment as PAML input file
%
% Syntax:  writepaml(aln,filename)
%
% Inputs:
%    aln         - Alignment structure
%    filename    - file name
%
% See also: READPHYLIP_I, WRITEPHYLIP_S

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if nargin < 2
    [filename, pathname, filterindex] = uiputfile( ...
       {'*.paml;*.pml', 'PAML Input Files (*.paml, *.pml)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];

	if (filterindex==1)
	if (isempty(find(filename=='.',1)))
		filename=[filename,'.pml'];
	end
	end
end

[n,m]=size(aln.seq);
p=1:n; q=1:m;
[NT,AA] = seqcode;

switch (aln.seqtype)
    case (3)	% Protein
	aln.seq(isnan(aln.seq))=i_getcode4gap('PROTEIN');
	Seq(p,q)=AA(aln.seq(p,q));
    otherwise	% nucleotides
	aln.seq(isnan(aln.seq))=i_getcode4gap('DNA');
	Seq(p,q)=NT(aln.seq(p,q));
end


fid = fopen(filename,'wt');
if fid == -1
   disp('Unable to open file.');
   return
end

fprintf(fid, ' %d %d\n',n,m);
mt = 1:60:size(Seq,2);
mt = cat(1,mt',size(Seq,2)+1);

for i=1:n
    fprintf(fid, '%s\n',aln.seqnames{i});
for j=1:length(mt)-1
    s=char(Seq(i,mt(j):mt(j+1)-1));
    fprintf(fid, '%s\n', s);
end
end
fclose(fid);

