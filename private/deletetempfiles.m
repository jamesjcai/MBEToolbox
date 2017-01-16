function deletetempfiles()
%DELETETEMPFILES

tempfiles={'2ML.dN' , '2ML.dS', '2ML.t', '2NG.dN', '2NG.dS', '2NG.t', 'seq_examples/CopyOfHCV_alignd.fas', ...
'seq_examples/CopyOfHCV_alignd_I.phy', 'seq_examples/CopyOfHCV_alignd_S.phy', 'infile', 'infile.dnd', 'outfile', ...
'outtree', 'rst', 'rst1', 'rub', 'seq_examples/unaligned_DNA.dnd', 'seq_examples/unaligned_Protein.dnd'};

for (k=1:length(tempfiles)),
	filename=char(tempfiles(k));
	if (exist(filename, 'file'))
	      delete(filename);
	end
end
