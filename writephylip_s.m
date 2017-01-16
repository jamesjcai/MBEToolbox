function writephylip_s(seq,filename)
%WRITEPHYLIP_S - Write alignment structure into a sequential PHYLIP formatted file
%
% Syntax:  writephylip_s(aln,'filename')
%
% Inputs:
%    aln         - Alignment structure
%    filename    - Sequencial PHYLIP file.
%
% See also: READPHYLIP_S, WRITEPHYLIP_I, READPHYLIP_I

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if isstruct(seq)
   aln=seq;
else
   aln.seq=seq;
end

if nargin < 2
    [filename, pathname, filterindex] = uiputfile( ...
       {'*.phylip;*.phy', 'Phylip Format Files (*.phylip, *.phy)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];

	if (filterindex==1)
		if (isempty(find(filename=='.')))
		filename=[filename,'.phy'];
		end
	end
end


[n,m]=size(aln.seq);
p=1:n; q=1:m;
[NT,AA] = seqcode;


if (isfield(aln,'seqtype'))
switch (aln.seqtype)
    case (3)	% Protein
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('PROTEIN');
	Seq(p,q)=AA(aln.seq(p,q));
    otherwise	% nucleotides
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('DNA');
	Seq(p,q)=NT(aln.seq(p,q));
end
else
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('DNA');
	Seq(p,q)=NT(aln.seq(p,q));
end

if (~isfield(aln,'seqnames'))
    for k=1:n
    aln.seqnames{k}=sprintf('Seq%d',k);
    end
end





fid = fopen(filename,'wt');
if fid == -1
   disp('Unable to open file.');
   return
end



fprintf(fid, [' %d %d\n'],n,m);
for i=1:n,
   name=i_name10(char(aln.seqnames(i)));
   fprintf(fid, ['%10s  %s\n'],name,char(Seq(i,:)));
end
fclose(fid);

