
trimAl 1.2rev57. Copyright (C) 2009. Salvador Capella-Gutierrez and Toni Gabaldón.

trimAl webpage: http://trimal.cgenomics.org

This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, the last available version.

Please cite:	Salvador Capella-Gutierrez, Jose M. Silla-Martinez and
            	Toni Gabaldon. trimAl: a tool for automated alignment 
            	trimming (2009).

Basic usage
	trimal -in <inputfile> -out <outputfile> -(other options).

Common options (for a complete list please see the User Guide or visit http://trimal.cgenomics.org):

    -h                       Print this information and show some examples.
    --version                Print the trimAl version.

    -in <inputfile>          Input file in several formats (clustal, fasta, NBRF/PIR, nexus, phylip3.2, phylip).
    -compareset <inputfile>  Input list of paths for the files containing the alignments to compare.
    -matrix <inpufile>       Input file for user-defined similarity matrix (default is Blosum62).

    -out <outputfile>        Output alignment in the same input format (default stdout). (default input format)
    -htmlout <outputfile>    Get a summary of trimal's work in an HTML file.

    -clustal                 Output file in CLUSTAL format
    -fasta                   Output file in FASTA format
    -nbrf                    Output file in NBRF/PIR format
    -nexus                   Output file in NEXUS format
    -mega                    Output file in MEGA format
    -phylip3.2               Output file in PHYLIP3.2 format
    -phylip                  Output file in PHYLIP/PHYLIP4 format

    -complementary           Get the complementary alignment.
    -colnumbering            Get the relationship between the columns in the old and new alignment.

    -select { n,l,m-k }      Selection of columns to be removed from the alignment. (see User Guide).
    -gt -gapthreshold <n>    1 - (fraction of sequences with a gap allowed).
    -st -simthreshold <n>    Minimum average similarity allowed.
    -ct -conthreshold <n>    Minimum consistency value allowed.
    -cons <n>                Minimum percentage of the positions in the original alignment to conserve.

    -nogaps                  Remove all positions with gaps in the alignment.
    -noallgaps               Remove columns composed only by gaps.

    -gappyout                Use automated selection on "gappyout" mode. This method only uses information based on gaps' distribution. (see User Guide).
    -strict                  Use automated selection on "strict" mode. (see User Guide).
    -strictplus              Use automated selection on "strictplus" mode. (see User Guide).
                             (Optimized for Neighbour Joining phylogenetic tree reconstruction).

    -automated1              Use a heuristic selection of the automatic method based on similarity statistics. (see User Guide).
                             (Optimized for Maximum Likelihood phylogenetic tree reconstruction).

    -resoverlap              Minimum overlap of a positions with other positions in the column to be considered a "good position". (see User Guide).
    -seqoverlap              Minimum percentage of "good positions" that a sequence must have in order to be conserved. (see User Guide).

    -w <n>                   (half) Window size, score of position i is the average of the window (i - n) to (i + n).
    -gw <n>                  (half) Window size only applies to statistics/methods based on Gaps.
    -sw <n>                  (half) Window size only applies to statistics/methods based on Similarity.
    -cw <n>                  (half) Window size only applies to statistics/methods based on Consistency.

    -sgc                     Print gap percentage count for columns in the input alignment.
    -sgt                     Print accumulated gap percentage count.
    -scc                     Print conservation values for columns in the input alignment.
    -sct                     Print accumulated conservation values count.
    -sfc                     Print compare values for columns in the selected alignment from compare files method.
    -sft                     Print accumulated compare values count for the selected alignment from compare files method.
    -sident                  Print identity statistics for all sequences in the alignemnt. (see User Guide).

