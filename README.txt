These scripts provide some basic functionality for BLAST-based analysis of plasmid sequences. For example output, and to cite these scripts, please use:
Stoesser et al (2015) bioRxiv doi: http://dx.doi.org/10.1101/030668

blast_plot.sh:
This takes a single reference plasmid sequence and uses BLASTn to compare it to a set of other sequences (typically de novo assemblies). Output is a figure (in pdf format) showing sequence similarity across the length of the plasmid reference. Default parameters produce a rainbow-coloured figure, which can be useful for distinguishing sequence identity values, but is not very pretty. Alternatively a single colour gradient can be specified. As an example, Figures 5 and 6 in the above cited paper were generated using the following parameters:
blast_plot.sh -l SAMPLES -r REFERENCE -f OUTFOLDER -o PDFFILE -R 90 -v darkred -S -d darkblue -C 0.2

reference_comparison.sh:
This takes a list of plasmid sequences and performs pairwise comparisons between them. Output is a matrix showing the percentage of the "reference" plasmid with similarity in the comparison plasmid. Note that even though each pair of plasmids is compared twice, the reference plasmid will be different in each case and therefore the matrix will not be symmetrical unless plasmids are the same length. 

Notes:
1. These scripts do not perform robust error checking. Please check inputs carefully.
2. Both scripts require the user to specify a folder for storing BLAST output files. The rationale for this is that with many sequences, blasting can take a long time. Therefore, results are stored such that e.g different colour parameters, or a modified sample list, can be specified without performing all reprocessing. In order to facilitate this, blast output file names are based on the base names of input files (i.e. without directory information). Therefore, DO NOT USE MULTIPLE INPUT FILES WITH THE SAME NAME. In addition, it has been noted that errors in specifying input correctly or processing interrupts can sometimes cause empty files to be produced, so please be aware of this and if in doubt, simply delete the contents of this folder to force reprocessing.
3. The additional .py and .R scripts are required for the above .sh scripts, but are not designed to be called directly.

Requirements:
BLAST+
python, Biopython
R (packages: gplots, argparse) - not required for reference_comparison.sh
