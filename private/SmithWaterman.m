% SmithWaterman:  Matlab implementation of the Smith-Waterman
%                 sequence alignment algorithm
%
% Scott F. Smith
% Department of Electrical and Computer Engineering
% Boise State University
% SFSmith@BoiseState.edu
%
% 20 February 2004

% Filename to read (fname.txt)
% Should contain two sequences in FASTA format
fname = 'test';

% Gap penalties
g = 5;
c = 3;

% List of allowed protein characters for PAM250
ProtChars = 'ABCDEFGHIKLMNPQRSTVWYZ';
LProtChars = 'abcdefghiklmnpqrstvwyz';

% Read FASTA format input file
PString1 = [];
PString2 = [];
fid = fopen([fname '.txt']);
% Throw away first (comment) line
FileLine = fgetl(fid);
% Read first sequence
while 1
  FileLine = fgetl(fid);
  if FileLine(1) == '>'
    break
  end
  PString1 = [PString1 FileLine];
end
while 1
  FileLine = fgetl(fid);
  if ~ischar(FileLine)
    break
  end
  PString2 = [PString2 FileLine];  
end
fclose(fid);

% List of allowed protein characters for PAM250
ProtChars = 'ABCDEFGHIKLMNPQRSTVWYZ';
LProtChars = 'abcdefghiklmnpqrstvwyz';

% The substitution matrix
PAM250 = [ ...
 2  0 -2  0  0 -4  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3  0
 0  2 -4  3  2 -5  0  1 -2  1 -3 -2  2 -1  1 -1  0  0 -2 -5 -3  2 
-2 -4 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0 -5 
 0  3 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4  3 
 0  2 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4  3  
-4 -5 -4 -6 -5  9 -5 -2  1 -5  2  0 -4 -5 -5 -4 -3 -3 -1  0  7 -5 
 1  0 -3  1  0 -5  5 -2 -3 -2 -4 -3  0 -1 -1 -3  1  0 -1 -7 -5 -1 
-1  1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0  2 
-1 -2 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1 -2 
-1  1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4  0 
-2 -3 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1 -3 
-1 -2 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2 -2 
 0  2 -4  2  1 -4  0  2 -2  1 -3 -2  2 -1  1  0  1  0 -2 -4 -2  1 
 1 -1 -3 -1 -1 -5 -1  0 -2 -1 -3 -2 -1  6  0  0  1  0 -1 -6 -5  0 
 0  1 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4  3  
-2 -1 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4  0 
 1  0  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3  0 
 1  0 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3 -1 
 0 -2 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2 -2 
-6 -5 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0 -6  
-3 -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10 -4 
 0  2 -5  3  3 -5 -1  2 -2  0 -3 -2  1  0  3  0  0 -1 -2 -6 -4  3 
];

SubMat = PAM250;

% Convert protein sequence ASCII characters to position in the substitution matrix
PIndex1 = zeros(size(PString1));
for i = 1:length(PString1)
  PIndex1(i) = [find(ProtChars == PString1(i)) find(LProtChars == PString1(i))];
end
PIndex2 = zeros(size(PString2));
for i = 1:length(PString2)
  PIndex2(i) = [find(ProtChars == PString2(i)) find(LProtChars == PString2(i))];
end

% Calculate scoring matrix
I = zeros(length(PString1)+1, length(PString2)+1);
D = zeros(length(PString1)+1, length(PString2)+1);
M = zeros(length(PString1)+1, length(PString2)+1);
trace = zeros(length(PString1)+1, length(PString2)+1);
for i = 2:length(PString1)+1
  for j = 2:length(PString2)+1
    I(i,j) = max(I(i-1,j)-c, M(i-1,j)-g);
    D(i,j) = max(D(i,j-1)-c, M(i,j-1)-g);
    M(i,j) = max([I(i-1,j-1) D(i-1,j-1) M(i-1,j-1)]);
    if I(i-1,j-1) == M(i,j)
      trace(i,j) = 1;
    end
    if D(i-1,j-1) == M(i,j)
      trace(i,j) = -1;
    end
    M(i,j) = M(i,j) + SubMat(PIndex1(i-1),PIndex2(j-1));
    M(i,j) = max(M(i,j), 0);    
  end
end
% Set traceback matrix on edges
trace(2:length(PString1)+1,2) = ones(length(PString1),1);
trace(2,2:length(PString2)+1) = -1 * ones(1,length(PString2));

% Traceback: 1 = up, 0 = diagonal, -1 = left
OutString1 = [];
OutString2 = [];
[vali maxi] = max(M(2:length(PString1)+1,length(PString2)+1));
[valj maxj] = max(M(length(PString1)+1,2:length(PString2)+1));
OutString1 = PString1(maxi);
OutString2 = PString2(maxj);
if maxi < length(PString1)
  OutString1 = [OutString1 PString1(maxi+1:length(PString1))];
  for i = maxi+1:length(PString1)
    OutString2 = [OutString2 '-'];
  end
end
if maxj < length(PString2)
  OutString2 = [OutString2 PString2(maxj+1:length(PString2))];
  for j = maxj+1:length(PString2)
    OutString1 = [OutString1 '-'];
  end
end
i = maxi + 1;
j = maxj + 1;
while i > 2 | j > 2
  if i > 2
    OutString1 = [PString1(i-2) OutString1];
  else
    OutString1 = ['-' OutString1];
  end
  if j > 2
    OutString2 = [PString2(j-2) OutString2];
  else
    OutString2 = ['-' OutString2]; 
  end
  if trace(i,j) > 0
    OutString2(1) = '-';
    nexti = i-1;
    nextj = j;
  end
  if trace(i,j) < 0
    OutString1(1) = '-';
    nexti = i;
    nextj = j-1;
  end
  if trace(i,j) == 0
    nexti = i-1;
    nextj = j-1;
  end
  i = nexti;
  j = nextj;
end

OutString1
OutString2

% Output aligned strings
fidOut = fopen([fname '_aligned.txt'],'w');
StartLoc = 0;
while StartLoc < length(OutString1)
  EndLoc = min([StartLoc+50 length(OutString1)]);
  fprintf(fid,'%.0f to %.0f\r\n',StartLoc+1,EndLoc);
  fprintf(fid,'%s\r\n',OutString1(StartLoc+1:EndLoc));
  fprintf(fid,'%s\r\n\r\n',OutString2(StartLoc+1:EndLoc));
  StartLoc = StartLoc + 50;
end
fclose(fidOut);