function [seq,header]=readfasta_largesingle(FASTAfilename)

seq=[];    
fidIn = fopen(FASTAfilename,'r');
header = fgetl(fidIn);
map = uint8([1 11 2 12 0 0 3 13 0 0 7 0 8 15 0 0 0 5 9 4 4 14 10 15 6 0 16]);
uint8a = uint8('a');
newLine = sprintf('\n');
blockSize = 2^20;
while ~feof(fidIn)
    charData = fread(fidIn,blockSize,'*char')';
    charData = strrep(charData,newLine,'');
    charData = lower(charData);
    intData = map(uint8(charData) - uint8a + 1);
    %fwrite(fidOut,intData,'uint8');
    seq=[seq, intData];
end
seq(seq>4)=5;
