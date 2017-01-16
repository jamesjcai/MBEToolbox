function [stats]=mapprun(aln)
%MAPPRUN - executes MAPP (MULTIVARIATE ANALYSIS OF PROTEIN POLYMORPHISM)
%
%mapprun(aln)
%[stats]=mapprun(aln)
%
%     stats.position - locations along a chromosome
%     stats.clr      - composite likelihood ratio
%     stats.alpha    - the strength of the sweep
%
% MAPP reads an alignment of protein sequences and a tree relating the 
% sequences. It then calculates the predicted impact of each potential SNP 
% at each position. The predictions are based on a set of scales (of 
% physicochemical properties) for which each amino acid has a numeric value.

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[exedir,dlgshown]=mbe_getprgmdir(sprintf('%s_prgmdir',mfilename));
if isempty(exedir)||dlgshown, return; end

oldpath=pwd;
cd(exedir);
[status]=snp_writesweepfinder(geno,mark,'input.tab');
if status~=1
    cd(oldpath);
    error('Error writing SweepFinder input file.');
end

% -f <path_and_alignment_file_name> Yes path to text file containing sequence alignment in Fasta format
% -t <path_and_tree_file_name> Yes path to text file containing tree with branch lengths
% -s <scales_to_be_used> No column numbers of scales to use, separated by colons (see above)
% -o <output_file_name> No path to output file

cmd=sprintf('java -jar MAPP.jar -f infile -t intree -o outfile.xls');
%fprintf('Running: %s\n\n',cmd);
system(cmd);


% -----------------

[vpos,vclr,valpha]=textread('output.tab','%f%f%f','headerlines',1);
if nargout<1
    vpos=vpos./1000000;
    figure;
    subplot(2,1,1)
    %plot(vpos,vclr,'LineSmoothing','on');
    plot(vpos,vclr);
    xlim([min(vpos(:)),max(vpos(:))]);
    ylabel('CLR (Composite Likelihood Ratio)')
    xlabel('Position (Mb)')
else
    stats=struct;
    stats.position=vpos;
    stats.CLR=vclr;
    stats.alpha=valpha;
end
cd(oldpath);

