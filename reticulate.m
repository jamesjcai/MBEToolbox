function [M] = reticulate(aln)
%RETICULATE - Calculats Compatibility Matrices
%A program for calculating and displaying compatibility matrices as an aid in
%determining reticulate evolution in molecular sequences.
%
%Compatibility Matrices
%In Jakobsen et.al. 1996 compatibility matrices were introduced for the
%study of reticulate evolution. The underlying idea is straighforward but
%nonetheless powerful. For every two parsimonious characters in a pre-given
%alignment, record in a dot-plot whether they are incompatible (i.e.
%conflicting) or compatible (not conflicting). Using two distinct colors
%(in our case black and white) for the two situations the recording is done
%in a symmetric matrix with the rows and columns successively labeled by
%the parsimonious sites of the alignment.
%
%REF:
%Jakobsen I.B., and Easteal S. (1996) A program for calculating and
%displaying compatibility matrices as an aid in determining reticulate
%evolution in molecular sequences, CABIOS, 12, 291-295.

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (isstruct(aln)), seq=aln.seq; else seq=aln; end

[seq,sites]=i_informative(seq);
[n,m]=size(seq);

if (m<2), error('RETICULATE:LessInformativeSites', 'Informative sites < 2'); end
disp('Display informative sites')
viewseq(seq);
M=zeros(m,m);


for i=1:m
    for j=i:m
        if i~=j
            su=unique(seq(:,[i,j]),'rows');
            M(i,j)=i_iscon(su);
            M(j,i)=M(i,j);
        end
    end
end

if (nargout<1)
	colormap([0 0 0; 1 1 1; .5 .5 .5]);
	X=M+1;
	X=X+eye(m)*2;
	image(X)
        axis square;
        % set('YTickLabel',num2str(sites));
end



function [y]=i_iscon(su)
	% Compatibility is calculated by comparing two informative sites. All distinct
	% combinations of characters at the two sites are found. For example, for the
	% five sequences below, the program had previously identified three informative
	%                             sites (they are repeated to the right of the
	%      2   6 8       2 6 8    alignment) then, when comparing sites 2 and 6,
	% A   ACGAACGTACGT   C C T    the pairs (C,C) (A,C) (A,T) are found, for
	% B   ACGTACGGACGT   C C G    2 and 8 they are (C,T) (C,G), (A,T) (A,G) and
	% C   AAGTACGTACGT   A C T    for 6 and 8, (C,T) (C,G) (T,G). Note that the
	% D   AAGTATGGACGC   A T G    order within the pairs is important. To find out
	% E   AAGTATGGACGT   A T G    if two sites are compatible, any pairs that have
	%                             a unique character in either position is
	% eliminated. So for sites 2 and 6 we can eliminate (C,C) since no other pair
	% has a C in the first position, and (A,T) since no other pair has a T in
	% the second position. The process is repeated and (A,C) can now be eliminated,
	% since it is the only pair left. Since all pairs are eliminated sites 2 and 6
	% are compatible. On the other hand, for sites 2 and 8, we cannot eliminate
	% any pair, as both C and A in the first position occurs in two pairs, and
	% similarly T and G in the second position. This means sites 2 and 8 are
	% incompatible. For sites 6 and 8, (C,T) and (T,G) can be eliminated - the
	% T is unique in each position. Thus 6 and 8 are also compatible.
	%
	%    2 6 8   The comparisons of sites are displayed in a matrix format.
	% 2  \ . M   The compatible comparisons are coloured white, the incompatible
	% 6  . \ .   are black. On the left, white squares are represented  by  .
	% 8  M . \   and  black by M . The diagonal (self with self) is \ .

	if (isempty(su)), y=1; end

	for n=1:2
	for k=1:4
	    x=find(su(:,n)==k);
	    if (length(x)==1), su(x(1),:)=[]; end
	end
	end
	if (isempty(su)), y=1; else y=0; end


function [S2,sites] = i_informative(seq)
	As=((sum(seq==1))>1);	% A
	Gs=((sum(seq==3))>1);	% G
	Cs=((sum(seq==2))>1);	% C
	Ts=((sum(seq==4))>1);	% T
	sites=find((As+Gs+Cs+Ts)>1);
	S2=seq(:,sites);
	[n,m]=size(S2);
	for k=m:-1:1
	      site=S2(:,k);
	      if length(unique(site))~=2
		    S2(:,k)=[];
	      end
	end
