function [s] = seqbyidt(transid,typeid)
%SEQBYIDT - only useful before Ensembl v50
%USAGE: seqbyidt('ENST00000307719',3)
%typeset={'genomic','cdna','coding','peptide','utr5','utr3'};

% Population Genetics & Evolution Toolbox, (C) 2007
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 2/23/2007

%transid='ENST00000307719';
%typeid=6;

typeset={'genomic','cdna','coding','peptide','utr5','utr3'};
typetxt=typeset{typeid};

% http://www.ensembl.org/info/website/archives/index.html

base='http://jul2008.archive.ensembl.org/Homo_sapiens/exportview?';
url=sprintf('type1=transcript&anchor1=%s&format=fasta&action=export&_format=Text&options=%s&options=&output=txt&submit=Continue+%3E%3E',transid,typetxt);
url=strcat(base, url);

	[s0,status0]=urlread(url);
	if ~(status0>0),
		disp('Unable to download data.');
		return;
	end
	[s] = strread(s0,'%s','delimiter','\n');
	s(1)=[];
	s=char(s);

