#!/bin/bash

scriptdir=$(dirname $0)

function usage {
    echo
    echo 'Compare (plasmid) references in a pairwise manner.'
    echo
    echo 'Usage:'
    echo 'bash reference_comparison.sh -l REFERENCES -f OUTFOLDER -o OUTFILE [-i IDTHRESH]'
    echo
    echo 'Arguments:'
    echo -e '\t-h, -u\tShow this help message and exit.'
    echo -e '\t-l\tFile containing list of tab separated reference name, fasta file pairs. Blast databases must have already been generated for each fasta file using default naming. Required.'
    echo -e '\t-f\tFolder for blast output files. Required.'
    echo -e '\t-o\tOutput file for comparison matrix. Required.'
    echo -e '\t-i\tMinimum sequence identity cutoff. Default: none (i.e. any blast hit is counted)'
}

while getopts "hul:r:f:o:i:" opt; do
    case "$opt" in
	h|u)
	    usage
	    exit 0
	    ;;
	\?)
	    argerror=true
	    ;;
	l)
	    reflist=$OPTARG
	    if [ ! -f $reflist ]; then
		echo "Error for -l, file does not exist: "$reflist
		argerror=true
	    fi
	    ;;
	f)
	    blastfolder=$(echo $OPTARG | sed 's/\([^/]\)$/\1\//g')
	    if [ ! -d $blastfolder ]; then
		echo "Error for -f, folder does not exist: "$blastfolder
		argerror=true
	    fi	    
	    ;;
	o)
	    outfile=$OPTARG
	    ;;
	i)
	    idthresh=$OPTARG
	    ;;
    esac
done


if [ ! $reflist ]; then
    echo "Missing required argument: -l"
    argerror=true
fi
if [ ! $blastfolder ]; then
    echo "Missing required argument: -f"
    argerror=true
fi
if [ ! $outfile ]; then
    echo "Missing required argument: -o"
    argerror=true
fi
if [ $argerror ]; then
    usage
    exit 1
fi


# Generate header line containing list of sequences used as references
(echo -ne reference'\t'
while read line; do
	refname=$(echo "$line" | cut -f1)
	echo -ne $refname'\t'
done < $reflist
echo) > $outfile  # if outfile already exists, this overwrites previous contents

# Loop over query sequences
while read qline; do
    qfasta=$(echo "$qline" | cut -f2)
    qname=$(echo "$qline" | cut -f1)
    echo -ne $qname'\t' >> $outfile  # Start each line with the name of the query sequence
  
    # create list of contig names and lengths, if it doesn't exist already
    contigsfile=$blastfolder$(basename $qfasta)-contigs.txt
    if [ ! -f $contigsfile ]; then
	python $scriptdir/fasta_length.py --fasta $qfasta --names > $contigsfile
    fi

    # Loop over reference sequences
    while read rline; do
	rfasta=$(echo "$rline" | cut -f2)
        rname=$(echo "$rline" | cut -f1)

	# do blast if output file doesn't already exist
	blastoutfile=$blastfolder$(basename $rfasta)-$(basename $qfasta).csv
	if [ ! -f $blastoutfile ]; then
	    blastn -query $qfasta -db $rfasta -outfmt 6 > $blastoutfile
	fi
	
	# calculate reference proportion for this comparison
	echo -ne $(python $scriptdir/parse_blast.py -i $blastoutfile -c $contigsfile -s --idthresh 1 | cut -f2)'\t' >> $outfile

    done < $reflist
    echo >> $outfile
done < $reflist

