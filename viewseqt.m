function viewseqt(seq)
%VIEWSEQT - View translate sequences
%
% Syntax:  viewseqt(seq)
%
% Inputs:
%    aln         - Alignment structure
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


[NT,AA] = seqcode;
[n,m]=size(seq);
if (n>1), seq=seq(1,:); n=1;
warning('Only works for the first sequence')
end


Seq=NT(seq);
    fprintf('\n');
	[ssiz,slen]=size(Seq);
	mt = 1:60:slen;
	mt = cat(1,mt',slen+1);


seqpx=AA(translateseq(seq));
Seqp=char(' '*ones(1,slen));
xk=1;
for (k=2:3:slen),
    Seqp(k)=seqpx(xk);
    xk=xk+1;
end


	for (j=1:length(mt)-1),
	for (i=1:n),
		if (j==1),
			fprintf(['%s\n'], lower(char(Seq(i,[mt(j):mt(j+1)-1]))));
			fprintf(['%s\n'], char(Seqp(i,[mt(j):mt(j+1)-1])));
        else
			fprintf(['%s\n'], lower(char(Seq(i,[mt(j):mt(j+1)-1]))));
			fprintf(['%s\n'], char(Seqp(i,[mt(j):mt(j+1)-1])));
		end
	end
		%fprintf('\n');
    end
