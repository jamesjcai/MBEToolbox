function writemega(aln,filename)
%WRITEMEGA - Write alignment structure into a MEGA formatted file
%
% Syntax:  writemega(aln,'filename')
%
% Inputs:
%    aln         - Alignment structure
%    filename    - FASTA formatted file (ASCII text file).
%
% See also: WRITEFASTA

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if nargin < 2
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.mega;*.meg', 'MEGA Format Files (*.mega, *.meg)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];
	if (filterindex==1)
		if (isempty(find(filename=='.')))
		filename=[filename,'.meg'];
		end
	end

end


[n,m]=size(aln.seq);
p=1:n; q=1:m;
[NT,AA] = seqcode;

switch (aln.seqtype)
    case (3)	% Protein
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('PROTEIN');
	% AA = 'ARNDCQEGHILKMFPSTWYV*-';
	% Seq(p,q)=AA(aln.seq(p,q));
	Seq = AA(aln.seq);
    otherwise	% nucleotides
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('DNA');
	% NT = ['ACGT-'];
	Seq(p,q)=NT(aln.seq(p,q));
end


fid = fopen(filename,'wt');
if fid == -1
   disp('Unable to open file.');
   return
end


mt = 1:60:size(Seq,2);
mt = cat(1,mt',size(Seq,2)+1);

fprintf(fid, '#Mega\n');
fprintf(fid, ' Title: MBEToolbox alignment output\n');
fprintf(fid, '\n');

for i=1:n,
   name=char(aln.seqnames(i));
   fprintf(fid, ['#%s\n'],name);

for (j=1:length(mt)-1),
   % fprintf(fid, ['%s\n'],char(Seq(i,:)));
   fprintf(fid, ['%s\n'],char(Seq(i,[mt(j):mt(j+1)-1])));
end
end

fclose(fid);