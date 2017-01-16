function [aln] = selectsitedlg(aln)
%SELECTSITEDLG - 

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (isstruct(aln)),
	seq=aln.seq;
else
	seq=aln;
end


items(1).name = 'Exclude';
items(1).default = 0;
items(1).exclusive = [2];
%items(1).indent = 1;

items(2).name = 'Include';
items(2).default = 1;
items(2).exclusive = [1];
%items(2).indent = 1;

items(3).name = 'Sites:';                 
items(3).default = 0;                       
items(3).values = {'23..25,'};        
items(3).help = 'Enter your site.';


title = 'Sequence type and genetic code';
% msg = sprintf(['Please select sequence type and genetic code']);
out = CSEFlagDialog(items, title);