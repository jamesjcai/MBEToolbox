function viewtree(tree)
%VIEWTREE - Displays a tree with njplot.

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

oldpath=pwd;
cdmbe;
cd 'addins';
cd 'njplot';
if (isunix), slas='./'; else slas=''; end


if (nargin<1),
	disp('Using an example tree')
	tree='((gorla-ECP:0.00664,(chimp-ECP:0.00578,((orang-ECP:0.02515,(((((human-EDN:0.00542,chimp-EDN:0.00312):0.00277,gorla-EDN:0.00365):0.01918,orang-EDN:0.02427):0.01979,macaq-EDN:0.07058):0.02561,tamar-EDN:0.11203):0.05274):0.02151,macaq-ECP:0.04586):0.00935):0.00064):0.00095,human-ECP:0.00095);'
end

%dirstr=pwd;
%if (ispc), 
% outtree = [dirstr,'\outtree'];
%else
% outtree = [dirstr,'/outtree'];
%end

fid = fopen('outtree','wt');
if fid == -1
   disp('Unable to open file.'); 
   return
end

fprintf(fid, ['%s'],tree);
fclose(fid);

if (ispc), 
	%system(['njplot.exe ', '"', outtree, '"']);
	system(['njplot.exe outtree']);
else
	system(['./njplot outtree']);
	%system(['./njplot ', '"', outtree, '"']);
end
cd(oldpath);
