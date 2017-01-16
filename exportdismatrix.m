function exportdismatrix(D,aln,header)
%EXPORTDISTANCEMATRIX - Export distance matrix to file

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if nargin<3
    header = [];
end;

if nargin<2
   nb=length(D);
   for ii = 1:nb
      names{ii} = ['Seq' int2str(ii)];
   end
   aln.seqnames=names;
end


ButtonName=questdlg('Do you want to save distance matrix?', ...
                    'Save as', ...
                    'Text','Excel','Word','Text');
switch ButtonName,
    case 'Text',
       writematrix(D,aln);
    case 'Excel',
      colnames = aln.seqnames;
      xlswrite(D,header,colnames);
    case 'Word',
      colnames = aln.seqnames;
      table2word(colnames,D);
end