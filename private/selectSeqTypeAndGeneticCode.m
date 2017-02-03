function [seqtype,geneticcode]=selectSeqTypeAndGeneticCode()

%prompt = sprintf('SEQUENCE TYPE:\n1 = Non-coding nucleotide\n2 = Coding nucleotide\n3 = Protein\n');
%in=inputdlg(prompt);
%seqtype=str2double(in{1});

%values = {'1 - Standard';'2 - Vertebrate Mithchondrial'; ...
%                   '3 - Yeast Mithchondrial';'4 - Mold Mithchondrial'; ...
%		   '5 - Invertebrate Mithchondrial';'6 - Protozoan Mithchondrial'; ...
%                   '7 - Coelenterate Mithchondrial';'8 - Mycoplasma'};

choice = questdlg('Select the type of input sequences:', ...
	'Select Sequence Type', ...
	'Noncoding Nucleotide','Coding Nucleotide','Protein','Noncoding Nucleotide');
% Handle response
switch choice
    case 'Noncoding Nucleotide'
        seqtype = 1;
    case 'Coding Nucleotide'
        seqtype = 2;
    case 'Protein'       
        seqtype = 3;
end
geneticcode=1;
