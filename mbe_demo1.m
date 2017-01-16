%% MBEToolbox DEMO - Alignment file IO
% Welcome to MBEToolbox.  This is a demonstration of
% MBEToolbox's file IO functions.
%
% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



cd(fileparts(which('mbe_demo1')))
% Working directory is changed.

%%
% First, let's just show how to read a FASTA formatted file
% into a MATLAB structure represents an alignment.

filename=fullfile('seq_examples','HCV_aligned.fas');
aln = readfasta(filename,2,1)

%%
% Now let's view the sequences in this alignment.
%
% VIEWSEQ displays sequences in PHYLIP format by default.

viewseq(aln)

%%
% If you use READFASTA function without specifying sequence type and/or
% genetic code, a dialogue will come up.

aln = readfasta(filename)

%%
% Make sure the alignment structure is a valid coding sequence alignment...

isvalidaln(aln,'CODING')

%%
% Of course, you can save the alignment into as a new FASTA formatted file.

filename=fullfile('seq_examples','CopyOfHCV_aligned.fas');
writefasta(aln,filename)

%%
% ...or you may save it as an interleaved PHYLIP formatted file,

filename1=fullfile('seq_examples','CopyOfHCV_aligned_I.phy');
writephylip_i(aln,filename1)

%%
% or as a sequential PHYLIP formatted file.

filename2=fullfile('seq_examples','CopyOfHCV_aligned_S.phy');
writephylip_s(aln,filename2)

%%
% Two functions are provided by MBEToolbox for read interleaved or
% sequential PHYLIP formatted file respectively.

aln2 = readphylip_i(filename1,2,1)
aln3 = readphylip_s(filename2,2,1)

%%
% If you have a set of unaligned DNA/Protein sequences in FASTA file,
% MBEToolbox, using its build-in CLUSTALW, can align these sequences
% and import the alignment.

AlnDNA = alignseqfile('DNA',['seq_examples',filesep,'unaligned_DNA.fas'])
AlnPro = alignseqfile('PROTEIN',['seq_examples',filesep,'unaligned_Protein.fas'])

%%
% MBEToolbox can align coding DNA sequences according to protein alignment.

AlnCDS = alignseqfile('CDS',['seq_examples',filesep,'unaligned_CDS.fas'],1)

%%
% Take a look of the alignment.

viewseq(AlnCDS)

%%
% Thank you for viewing this introduction to MBEToolbox file IO funcation.