function [good,aln,filename]=checkcdsfile(filename,seqtype,geneticcode)
%CHECKCDSFILE - 

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if nargin < 1
    [filename, pathname] = uigetfile( ...
       {'*.fasta;*.fas', 'FASTA Format Files (*.fasta, *.fas)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a FASTA file');
	if ~(filename), aln=[]; good=0; return; end
	filename=[pathname,filename];
end

if nargin < 2
	seqtype=2; geneticcode=1;
end
good=1;


[aln] = readfasta(filename,seqtype,geneticcode);
if (isempty(aln)),
      good=0;
      return;
end


[n,m]=size(aln.seq);

x=[];

disp('Checking CDS files')

for (k=1:n),
      seq=aln.seq(k,:);
      seq(seq>4)=[];
      seqlen=length(seq);
      seqname=char(aln.seqnames(k));
     if (mod(seqlen,3)>0), 
        fprintf('Length of %d bp cannot divide by 3! -- %s\n',seqlen, seqname); 
	x=[x,k];
	good=0;
     end
end
      aln.seq(x,:)=[];
      aln.seqnames(x)=[];
