function runmapp(aln,tree)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if nargin<2
   [tree] = gettreedlg(aln);
end

oldpath=pwd;
cdmbe;
cd 'addins';
cd 'MAPP';
if (isunix), slas='./'; else slas=''; end

writefasta(aln,'infile');
fid=fopen('intree','w');
fprintf(fid,'%s\n',tree);
fclose(fid);
%java -jar MAPP.jar -f LacI_Alignment.fa -t LacI.tree -o LacI_output.xls
system(['java -jar MAPP.jar -f infile -t intree -o outfile']);

%[Position	Column Score	Column p-value	Alignment	Gap Weight	Over Gap Weight Threshold	Hydropathy	Polarity	Charge	Volume	Free energy alpha	Free energy beta	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	Good Amino Acids	Bad Amino Acids
%    names,types,x,y,answer] = textread('outfile','%s%s%f%d%s','headerlines',1);

cd(oldpath);


function data=dataimport
% Import data from outfile
% Automatically generated 22-Apr-2009

% Define parameters
fileName='C:\myprojects\mbetoolbox\addins\MAPP\outfile';
numHeaderLines=1;
formatString='%*q%f%f%q%d8%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q';
numRows=100;

% Read data from file
fid=fopen(fileName,'rt');
data=textscan(fid,formatString,numRows,'headerlines',numHeaderLines,'delimiter','\t');
status=fclose(fid);

