function frametrans(seq)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

disp('5''3'' Frame 1')
viewseqt(seq(1,:))
fprintf('\n');

disp('5''3'' Frame 2')
viewseqt(seq(1,2:end-2))
fprintf('\n');

disp('5''3'' Frame 3')
viewseqt(seq(1,3:end-1))
fprintf('\n');

