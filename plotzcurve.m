function plotzcurve(aln,action)
%PLOTZCURVE - Z curve plotter
%
%REFERENCE:
% C.T. Zhang & R. Zhang, (1991) Analysis of distribution of bases in the
%      coding sequences by a diagrammatic technique.
%      Nucl. Acids Res., 19, 6313-6317.
%
% R. Zhang & C.T. Zhang, (1994) Z Curves, an Intuitive Tool for Visualizing
%      and Analyzing DNA sequences.
%      J. Biomol. Struc. Dynamics 11, 767-782.
%
% Syntax: plot_z_curve(aln,'action');
%
% Inputs:
%    aln      - alignement/sequence object, i.e., a matrix/vector representation of
%               DNA sequence
%    action   - (optional) '3d'|'gc'|'3frame'
%
%
% Other m-files required: ZCURVE.m
% Subfunctions: none
% MAT-files required: none
%
% See also: ZCURVE

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if nargin<2,
    action='3d';
end;


if strcmp(action,'3d'),
	Seq=aln.seq(1,:);
        [X,Y,Z,Z1,AT,GC]=zcurve(Seq);
	plot3(X,Y,Z,'-');
	% plot3(M(1,:),M(2,:),M(3,:),'-','LineWidth',1);
	%hold on
	%fnplt(cscvn(M),'r',2)
	%hold off
	% set(gca,'XTick',[],'YTick',[],'ZTick',[])
	ylabel('Y'); 	xlabel('X');  	zlabel('Z');
	axis vis3d
	box on

elseif strcmp(action,'gc'),
        Seq=aln.seq(1,:);
	[X,Y,Z,Z1,AT,GC]=zcurve(Seq);
	subplot(2,1,1),plot(Z)
	title('Z curve (A+T)-(C+G)');
	xlabel('n (bp)'); ylabel('Z (bp)');

	subplot(2,1,2),plot(Z1)
	title('Cumulative GC profile (Z'' curve)');
	xlabel('n (bp)'); ylabel('Z'' (bp)');


elseif strcmp(action,'3frames'),


	Aln1=extractpos(aln,1); Seq1=Aln1.seq(1,:); [X1,Y1,Z1,Z11,AT1,GC1]=zcurve(Seq1);
	aln2=extractpos(aln,2); Seq2=aln2.seq(1,:); [X2,Y2,Z2,Z12,AT2,GC2]=zcurve(Seq2);
	Aln3=extractpos(aln,3); Seq3=Aln3.seq(1,:); [X3,Y3,Z3,Z13,AT3,GC3]=zcurve(Seq3);
	subplot(3,1,1),plot(Z1)
	title('Cumulative GC profile (Z'' curve) of frame 1');
	subplot(3,1,2),plot(Z2)
	title('Cumulative GC profile (Z'' curve) of frame 2');
	ylabel('Z'' (bp)');
	subplot(3,1,3),plot(Z3)
	title('Cumulative GC profile (Z'' curve) of frame 3');
    	xlabel('n (bp)');
elseif strcmp(action,'info');
    helpwin(mfilename);
end;    % if strcmp(action, ...