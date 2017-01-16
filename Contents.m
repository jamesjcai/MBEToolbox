% MBEToolbox -- Molecular Biology & Evolution Toolbox (old files)
%
% Data type and file IO
%  codonise64                - Codonises sequence(s)
%  copyalnheader             - Copy header information of an alignment structure into another
%  encodeseq                 - Converts nucleotide from letters to integer
%  seqcode                   - Return matrix for mapping sequence letters to integers
%  codontable                - Return codon matrix and genetic tables matrix
%  readfasta                 - Reads data from a FASTA formatted file into a MATLAB structure
%  readphylip_i              - Reads data from an interleaved PHYLIP formatted file
%  readphylip_s              - Reads data from a sequential PHYLIP formatted file
%  printmatrix               - Prints to the screen a printout of the matrix
%  viewseq                   - View sequences in alignment
%  writefasta                - Write alignment structure into a FASTA formatted file
%  writephylip_i             - Write alignment structure into an interleaved PHYLIP formatted file
%  writephylip_s             - Write alignment structure into a sequential PHYLIP formatted file
%  writematrix               - Writes data in tabular form to the file system.
%  writeshadedhtml           - Write shaded alignment into HTML file
%
% Basic statistic and sequence manipulation
%  countntchange             - Count nucleotide changes in two DNA sequences
%  countaachange             - Count amino-acid changes in two protein sequences
%  countseqppq               - Counts transitions (P1 and P2) and transversion (Q)
%  countseqpq                - Counts transition (P) and transversion (Q) for the given sequence pair S
%  gc4                       - Counts GC content at fourfold degenerate sites
%  cai                       - Calculates the Codon Adaptation Index (CAI)
%  rscu                      - Calculates the Relative Synonymous Codon Usage (RSCU) value
%  codonusage                - Counts codon usage
%  codonvolatility           - Calculates codon volatility
%  countdegeneratesites      - Counts degenerate sites in two aligned DNA sequences.
%  countinvariablesites      - Counts invariable sites.
%  ntcomposition             - Counts nucleotide composition
%  aacomposition             - Counts AA composition
%  countsegregatingsites     - Count segregating sites
%  extractdegeneratesites    - Extract 0-, 2- and 4-fold degenerate sites
%  extractinformativesites   - Extract informative sites
%  extractinvariablesites    - Extract invariable sites
%  extractpos                - Extract coding position 1, 2, 3 or 1 and 2
%  extractsegregatingsites   - Extract segregating sites
%  getsynnonsyndiff          - Return matrices of Syn- Nonsyn- differences between codons
%  getsynnonsynsites         - Return matrices of Syn- Nonsyn- sites of codons
%  hasgap                    - Check if alignment contains gap
%  rmcodongaps               - Remove codons with gaps 
%  rmgaps                    - Remove gaps in alignment
%  revcomseq                 - Return reverse complement of sequences
%  revseq                    - Return reverse strand of nucleotide sequences
%  translatealn              - Translate coding DNA sequence into protein sequence in an alignment
%  translateseq              - Translate coding DNA sequence into protein.
%  smithwaterman             - Local alignment by Smith & Waterman algorithm
%  needlemanwunsch           - Gobal alignment of two proteins.
%  zscoreprotaln             - Z score of protein alignment
%  alignseqfile              - Align sequence file
%  clustalw                  - Mulitple sequence alignment
%  karlinsig                 - Return the Karlin genomic signatures for a given sequence
%  karlinsigdiff             - Return the difference of Karlin genomic signatures between two given sequences
%  buildpssm                 - Builds a position specific scoring matrix (PSSM)
%  compcomp                  - Calculates the compositional complexity of sequence
%  cataln                    - Concatenate alignments
%  randcut                   - Cuts genome DNA randomly into shot fragment
%
% Genetic distances
%  dc_ng86                   - Compute syn. non-syn. substitutions rates by Nei-Gojobori method
%  dc_li85                   - Compute syn. non-syn. substitutions rates by Li85 method
%  dc_li93                   - Computes Syn. non-syn. substitutions rates by Li93 method
%  dc_ml                     - Estimate dS and dN by using maximum likelihood algorithm
%  dn_ntdiff                 - Number of different nucleotides and gaps between two sequences
%  dn_pdist                  - p-distances (nucleotide)
%  dn_ntfreqdist             - Euclidean distances between nucleotide frequencies
%  dn_jc                     - Jukes-Cantor Distance
%  dn_k2p                    - Kimura 80 (2-parameter) Distance
%  dn_tajima_nei84           - Tajima & Nei 84 Distance
%  dn_tamura92               - Tamura 92 Distance
%  dn_f84                    - Felsenstain 84 Distance
%  dn_hky                    - Hasegawa, Kishino and Yano 85 (HKY) Distance
%  dn_gtr                    - Distance of GTR model
%  dn_jin_nei90              - Jin-Nei Gamma distance
%  dn_logdet                 - Log-det (paralinear) distance
%  dn_ntfreqdist             - Euclidean distances between nucleotide frequencies 
%  dp_aadiff                 - Uncorrected distance of protein sequences
%  dp_pdist                  - p-distances (amino acid)
%  dp_poisson                - Poisson Correction (PC) distance (amino acid)
%  dp_gamma                  - Gamma distance (amino acids)
%  dp_dayhoff                - Dayhoff Distance
%  dp_jtt                    - JTT Distance
%  dp_wag                    - WAG Distance
%  estimatefreq              - Estimates base frequencies of the given sequence(s)
%  estimatekappa             - Estimates kappa for given sequence pair
%  estimateyn00kappa         - Estimate transition/transversion rate ratio (kappa) by method of YN00
%  gammadistrib              - Discrete Gamma model of rate heterogeneity
%  invdistrib                - Invariable sites model of rate heterogeneity
%
% Phylogeny and likelihood
%  mbe_dnaml                 - M-file for mbe_dnaml.fig
%  mbe_dnapars               - M-file for mbe_dnapars.fig
%  mbe_proml                 - M-file for mbe_proml.fig
%  mbe_protpars              - M-file for mbe_protpars.fig
%  plotupgma                 - Performs UPGMA on distance matrix and produces plot of dendrogram
%  plotnjtree                - Performs neighbor joining (NJ) on distance data and plot NJ tree
%  runaddin                  - Run add-in Phylip commands
%  likelidist                - Estimates log likelihood of branch length (distance)
%  seqpairlikeli             - Estimates log-likelihood of branch length (distance)
%  likelitree                - Estimates log likelihood of a tree
%  model_nt                  - Returns a structure of nucleotide substitution model
%  modeljc                   - Returns a structure of model JC
%  modelk2p                  - Returns a structure of model K2P
%  modelf81                  - Returns a structure of model F81
%  modelhky                  - Returns a structure of model HKY
%  modelgtr                  - Returns a structure of model GTR
%  model_aa                  - Returns a structure of amino-acid substitution model
%  modeldayhoff              - Returns a structure of model Dayhoff
%  modeljtt                  - Returns a structure of model JTT
%  modelwag                  - Returns a structure of model WAG
%  randtre                   - Generates a random tree with n OTUs
%  readnewick                - Reads tree from standart file in Newick format
%  sitepattern              - Infers site patterns for a given alignment
%  optimlikelidist           - Optimises distance of given sequence pair under a model
%  optimlikelidistk2p        - Optimises distance and kappa under a K2P model
%  optimseqpairlikeli        - Optimises distance of given sequence pair under a model
%  parsetree                 - Parses string of tree in Newick format
%  composeQ                  - Computes normalized rate matrix Q from R matrix (general reversible model)
%  compequtest               - Compositional equilibrium test
%
% Graph and plot
%  plotcorresp               - Performs correspondence analysis and plot result
%  plotdiplomo               - Plot DiPloMo graph
%  plotdistvstrans           - Plot genetic distances (NG86) vs. transitions and transversions
%  plotkavsks                - Plot Ka vs. Ks
%  plotntcomposition         - Plot nucleotide composition
%  plotslidingwin            - Performs sliding window analysis on a nucleotide sequence
%  plotslidingwinkaks        - Plot cumulative Ka and Ks curves
%  plotzcurve                - Z curve plotter
%  zcurve                    - Return Z curve components
%
% Polymorphism
%  reportpolysites           - Report polymorphic sites
%  tajima89test              - Tajima's Test of Neutrality

% MBEGUI                      - MBEToolbox GUI
% MBEDEMO                     - DEMO1: Alignment file IO 
%                             - DEMO2: Basic sequence statistics
%                             - DEMO3: Genetic distances estimation
%                             - DEMO4: Phylogenetic inferences
%
%% $Date: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $