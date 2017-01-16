function writeshadedhtml(aln,filename,level)
%WRITESHADEDHTML - Write shaded alignment into HTML file
%This funciton produces a display that emphasizes the degree of conservation
%in each column in the alignment. This focuses attention on which parts of the
%sequences are least tolerant of change, or where a change in the sequence is
%most likely to change the structure or function of the molecule represented
%by the sequence. Up to four different Shading Levels, LEVEL, of conservation
%can be distinguished.  The degree of conservation at each level can be set
%by the user in the Configuration dialog.  Default levels are set to complete
%or 100 percent conserved, 80 percent or greater conserved, 60 percent or
%greater conserved, and less than 60 percent conserved.
%
%This menu item toggles the ‘Enable Similarity Groups’ flag, also found in
%the Configuration dialog on the Similarity Groups page and as a button on the tool bar. If this flag is set,
%a Similarity Group will be used if at least any two of the characters found in the Similarity Group are
%found in the alignment column.  A Similarity Group is a set of sequence residues that are to be treated as if they are equivalent.  One consistent basis for selecting the members of a Similarity Group is to select a set of residues where the score for any pair of residues in the set is greater than zero in the current scoring table.  The default Similarity Groups for each scoring table have been selected this way. Thus, these Similarity Groups represent what are usually called conservative substitutions, evolutionary substitutions that are evidence in favor of the sequences being homologous (related through a common evolutionary ancestor).
%
% Syntax:  writeshadedaln(aln,'filename')
%
% Inputs:
%    aln         - Alignment structure
%    filename    - Name of shaded alignment file.
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


if (nargin<3), level=2; end

if nargin < 2
    [filename, pathname, filterindex] = uiputfile( ...
       {'*.html;*.htm', 'HTLM Format Files (*.html, *.htm)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];

	if (filterindex==1)
	if (isempty(find(filename=='.')))
		filename=[filename,'.htm'];
	end
	end
end

[n,m]=size(aln.seq);
p=1:n; q=1:m;
[NT,AA] = seqcode;

switch (aln.seqtype)
    case (3)	% Protein
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('PROTEIN');
	Seq(p,q)=AA(aln.seq(p,q));
	seqdeg=i_degseq(aln.seq,level);
	seqrel=aln.seq;
    otherwise	% nucleotides
	aln.seq(find(isnan(aln.seq)))=i_getcode4gap('DNA');
	Seq(p,q)=NT(aln.seq(p,q));
	seqdeg=aln.seq;
	seqrel=aln.seq;
end


pos1=i_primarylevel(seqdeg);
[pos2,exceptpos2]=i_secondarylevel(seqdeg);
[pos3,exceptpos3]=i_tertiarylevel(seqdeg);




fid = fopen(filename,'wt');
if fid == -1
   disp('Unable to open file.');
   return
end
fprintf(fid,'<pre>\n');
fprintf(fid, [' %d %d\n'],n,m);
mt = 1:60:size(Seq,2);
mt = cat(1,mt',size(Seq,2)+1);
star='*';

for (j=1:length(mt)-1),
	fprintf(fid,'%10s   %10s%10s%10s%10s%10s%10s\n', ' ', '*', num2str(20*((j-1)*3+1)), '*', num2str(20*((j-1)*3+2)), '*', num2str(20*((j-1)*3+3)));
for (i=1:n),

	t=char(Seq(i,[mt(j):mt(j+1)-1]));  s='';
	for (k=1:length(t)),
	     pos=(mt(j)-1)+k;
	     if (ismember(pos,pos1))
		     s=[s, sprintf('<span style="background: black; color: white">%s',t(k)),'</span>'];
             elseif (ismember(pos,pos2))
	             if (exceptpos2(i,pos)==1)
                        s=[s, t(k)];
		     else
		        s=[s, sprintf('<span style="background: gray; color: white">%s',t(k)),'</span>'];
		     end
             elseif (ismember(pos,pos3))
    	             if (exceptpos3(i,pos)==1)
                        s=[s, t(k)];
		     else
		        s=[s, sprintf('<span style="background: silver">%s',t(k)),'</span>'];
		     end
             else
		     s=[s, t(k)];
	     end
	end


	if (j==1)
		name=char(aln.seqnames(i));
		fprintf(fid, ['%10s : %s :%5d\n'], i_name10(name), s, min(j*60,size(Seq,2)));
	else
		name=char(aln.seqnames(i));
		fprintf(fid, ['%10s : %s :%5d\n'], i_name10(name), s, min(j*60,size(Seq,2)));
		% fprintf(fid, ['%10s %s\n'],' ', s);
	end
end
	fprintf(fid,'\n');
end
fprintf(fid,'</pre>\n');
fclose(fid);

%web(['file://', filename],'-browser')



%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%


function [pos] = i_primarylevel(s)
[n,m]=size(s);
pos=zeros(1,m);
for j=1:m,
	minnt=min(s(:,j));
	maxnt=max(s(:,j));
	pos(1,j)=(minnt==maxnt);
end
pos=find(pos==1);


function [pos,exceptpos] = i_secondarylevel(s)
[n,m]=size(s);
%exceptpos=zeros(n,m);
exceptpos=(s>20);
pos=zeros(1,m);
for j=1:m,
	minnt=min(s(:,j));
	maxnt=max(s(:,j));
	x=sum(s(:,j)==minnt);
	y=sum(s(:,j)==maxnt);
	if (x<y)
		exceptpos(:,j)=(s(:,j)==minnt);
	else
		exceptpos(:,j)=(s(:,j)==maxnt);
	end


        if (sum([x y])==n)
	     %pause
	     xx=max(x,y)./(max(x,y)+min(x,y));
	     if (xx>=0.8),
	         pos(1,j)=1;
	     end
        end
end
pos=find(pos==1);



function [pos,exceptpos] = i_tertiarylevel(s)
[n,m]=size(s);
%exceptpos=zeros(n,m);
exceptpos=(s>20);
pos=zeros(1,m);
for j=1:m,
	minnt=min(s(:,j));
	maxnt=max(s(:,j));
	x=sum(s(:,j)==minnt);
	y=sum(s(:,j)==maxnt);
	if (x<y)
		exceptpos(:,j)=(s(:,j)==minnt);
	else
		exceptpos(:,j)=(s(:,j)==maxnt);
	end


        if (sum([x y])==n)
	     %pause
	     xx=max(x,y)./(max(x,y)+min(x,y));
	     if (xx>=0.6),
	         pos(1,j)=1;
	     end
        end
end
pos=find(pos==1);




function [seqdeg] = i_degseq(seqrel,deglevel)
%I_DEGSEQ -

seqdeg=seqrel;
if (nargin<2), deglevel=2; end
switch (deglevel)
    case (1)
	setaa={encodeseq('DEHKRNQST',3), encodeseq('LIVMFYWAGCP',3)};
    case (2)
	setaa={encodeseq('DEHKR',3), encodeseq('NQST',3), encodeseq('LIVMFYW',3),encodeseq('AG',3), encodeseq('CP',3)};
    case(3)
	setaa={encodeseq('DE',3), encodeseq('HKR',3), encodeseq('NQ',3),encodeseq('ST',3), encodeseq('LIV',3),...
	encodeseq('FYW',3), encodeseq('AG',3), encodeseq('MC',3), encodeseq('P',3)};
end
for (k=1:length(setaa)),
      x=setaa{k};
      x1=x(1);
      seqdeg(ismember(seqdeg,x))=x1;
end


% The Program Property Defaults
% Level 1
%               polar                              non-polar
%         D,E,H,K,R,N,Q,S,T                   L,I,V,M,F,Y,W,A,G,C,P


% Classification by charge was made by dividing the
% amino acids into three categories: positive (R, H, K),
% negative (D, E), and uncharged (A, N, C, Q, G, I, L, M,
% F, P, S, T, W, Y, V).


% Level 2
%       polar          polar       (non-polar)                      not
%      charged        uncharged    hydrophobic       small        grouped
%     D,E,H,K,R       N,Q,S,T       L,I,V,M,F,Y,W      A,G           C,P
%
% Level 3
%  negatively   positively                                                                          other
%  charged      charged    amide  alcohol  aliphatic  aromatic  small  sulfur
%   D,E          H,K,R          N,Q      S,T      L,I,V        F,Y,W     A,G    M,C    P




% For nucleic acid sequences:
% Grouping 1 = default
%           Purine = Large            Pyrimidine = Small
%               A,G                               C,T or U
% Grouping 2
%           Strong Hydrogen Bonding          Weak Hydrogen Bonding
%                    G,C                                         A,T or U
% Grouping 3

% A, C, G, TU, one in each group. When combined with the property mode style
% setting of Shade All, this group can be used to highlight all the difference
% residue types of a DNA project.  This is useful in making 'islands' of different
% types of nucleotide patterns stand out from regions that are more random in
% composition.
%
% Classification by charge was made by dividing the
% amino acids into three categories: positive (R, H, K),
% negative (D, E), and uncharged (A, N, C, Q, G, I, L, M,
% F, P, S, T, W, Y, V).
%
% Classification by volume and polarity was made by
% dividing the amino acids into six categories: special (C),
% neutral and small (A, G, P, S, T), polar and relatively
% small (N, D, Q, E), polar and relatively large (R, H, K),
% nonpolar and relatively small (I, L, M, V), and nonpolar
% and relatively large (F, W, Y).