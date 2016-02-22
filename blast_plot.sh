#!/bin/bash

scriptdir=$(dirname $0)

function usage {
    echo
    echo 'Generate a heatmap based on blast identity'
    echo
    echo 'Usage:'
    echo 'bash blast_plot.sh -l SAMPLES -r REFERENCE -f OUTFOLDER -o PDFFILE [-c -S -g -m TEMPFILE -b BINSIZE -n NEWSAMPLEORDER -C CHARACTERSIZE -P PDFSIZE -R RAINBOWSTART -p PRESENCEFILE -k SCALEINCREMENT -v CONTIGSEPCOLOUR -d SINGLECOLOUR]'
    echo
    echo 'Arguments:'
    echo -e '\t-h, -u\tShow this help message and exit.'
    echo -e '\t-l\tFile containing list of tab separated sample name, blast database pairs. Required.'
    echo -e '\t-r\tReference fasta file. Required.'
    echo -e '\t-f\tFolder for blast output files. Required.'
    echo -e '\t-m\tTemporary file name for plotting matrix. Default: tmp.csv'
    echo -e '\t-o\tOutput file name (pdf format). Required.'
    echo -e '\t-b\tBin size. Must be a positive integer. Default: 100'
    echo -e '\t-c\tReorder samples according to automatic clustering. Default is to output in the same order as in the input file'
    echo -e '\t-n\tName of output file for writing sample order. Only useful in practice if reordering is performed.'
    echo -e '\t-C\tCharacter size expansion for sample names.'
    echo -e '\t-P\tPdf dimension expansion factor.'
    echo -e '\t-R\tPercentage identity value at which rainbow colours should start. Default: 80%'
    echo -e '\t-p\tIf specified, name of file to write presence / absence values to. Presence defined as >=99% sequence identity over >=80% of reference length.'
    echo -e '\t-k\tScale bar increments, in kb. Default: 10'
    echo -e '\t-v\tColour to use for vertical line separating contigs (red, black etc). Default: black'
    echo -e '\t-S\tFor multi-fasta references, use a single scale instead of restarting at each contig break.'
    echo -e '\t-d\tUse a single specified colour (gradient) instead of a rainbow colour scheme (blue, black etc).'
    echo -e '\t-g\tDo not use a gradient, but rather plot presence/absence according to given threshold (RAINBOWSTART). Must specify colour to use (SINGLECOLOUR).'
}

matrixfile="tmp.csv"
binsize=100
blastlengthcutoff=0

while getopts "hul:r:f:m:o:b:cn:C:P:R:p:k:v:Sd:g" opt; do
    case "$opt" in
	h|u)
	    usage
	    exit 0
	    ;;
	\?)
	    argerror=true
	    ;;
	l)
	    velvetlist=$OPTARG
	    if [ ! -f $velvetlist ]; then
		echo "Error for -l, file does not exist: "$velvetlist
		argerror=true
	    fi
	    ;;
	r)
	    ref=$OPTARG
	    if [ ! -f $ref ]; then
		echo "Error for -r, file does not exist: "$ref
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
	m)
	    matrixfile=$OPTARG
	    ;;
	o)
	    outpdf=$(echo "$OPTARG".pdf | sed 's/.pdf.pdf$/.pdf/g') # Append .pdf if not provided
	    ;;
	b)
	    binsize=$OPTARG
	    ;;
	c)
	    reorder=true
	    ;;
	n)
	    orderfile=$OPTARG
	    ;;
	C)
	    cexRow=$OPTARG
	    ;;
	P)
	    pdfExp=$OPTARG
	    ;;
	R)
	    rainstart=$OPTARG
	    ;;
	p)
	    presencefile=$OPTARG
	    ;;
	k)
	    scaleinc=$OPTARG
	    ;;
	v)
	    contigcol=$OPTARG
	    ;;
	S)
	    singlescale=true
	    ;;
	d)
	    singlecol=$OPTARG
	    ;;
	g)
	    nogradient=true
	    ;;
    esac
done

if [ ! $velvetlist ]; then
    echo "Missing required argument: -l"
    argerror=true
fi
if [ ! $ref ]; then
    echo "Missing required argument: -r"
    argerror=true
fi
if [ ! $blastfolder ]; then
    echo "Missing required argument: -f"
    argerror=true
fi
if [ ! $outpdf ]; then
    echo "Missing required argument: -o"
    argerror=true
fi
if [ $argerror ]; then
    usage
    exit 1
fi

	
# create list of contig names and lengths, if it doesn't exist already
contigsfile=$blastfolder$(basename $ref)-contigs.txt
if [ ! -f $contigsfile ]; then
    python $scriptdir/fasta_length.py --fasta $ref --names > $contigsfile # $ref must be valid fasta file
fi

echo -n sample > $matrixfile
cut -f2 $contigsfile | while read contiglen; do
    echo -n ,
    bins=$(($contiglen / $binsize))
    eval echo -n {1..$bins} | sed 's/ /,/g'
done >> $matrixfile
echo >> $matrixfile

# if presencefile is specified, and it already exists, delete it
if [ $presencefile ]; then
    if [ -f $presencefile ]; then
	rm $presencefile
    fi
fi

# for each sample, perform blast, generate binned values for plotting, and optionally determine presence / absence
while read line; do
    samplename=$(echo "$line" | cut -f1)
    blastdb=$(echo "$line" | cut -f2) # Blast db must already exist
    blastoutfile=$blastfolder$(basename $blastdb)-$(basename $ref).csv
    # if blast output file already exists, use this instead of reprocessing
    if [ ! -f $blastoutfile ]; then
	blastn -db $blastdb -query $ref -outfmt 6 > $blastoutfile
    fi
    binnedfile=$blastoutfile-$binsize.csv

    # Parse blast output to generate binned file for plotting and determine presence values, if needed. If binned file already exists, use this instead of reprocessing.
    if [ ! -f $binnedfile ]; then
	if [ $presencefile ]; then
	    (echo -ne $samplename'\t'; python $scriptdir/parse_blast.py -i $blastoutfile -c $contigsfile -b $binsize -o $binnedfile -s | cut -f3) >> $presencefile
	else
	    python $scriptdir/parse_blast.py -i $blastoutfile -c $contigsfile -b $binsize -o $binnedfile
	fi
    elif [ $presencefile ]; then
	(echo -ne $samplename'\t'; python $scriptdir/parse_blast.py -i $blastoutfile -c $contigsfile -s | cut -f3) >> $presencefile
    fi
    (echo -n $samplename,
	cat $binnedfile) >> $matrixfile

done < $velvetlist

# call R script for plotting, with appropriate arguments
rargs="--matfile $matrixfile --pdffile $outpdf --scale $binsize --contigs $contigsfile"
if [ $reorder ]; then
    rargs="$rargs --reorder"
fi
if [ $orderfile ]; then
    rargs="$rargs --orderfile $orderfile"
fi
if [ $cexRow ]; then
    rargs="$rargs --cexRow $cexRow"
fi
if [ $pdfExp ]; then
    rargs="$rargs --pdfExp $pdfExp"
fi
if [ $rainstart ]; then
    rargs="$rargs --rainbowstart $rainstart"
fi
if [ $scaleinc ]; then
    rargs="$rargs --scaleinc $scaleinc"
fi
if [ $contigcol ]; then
    rargs="$rargs --contigcol $contigcol"
fi
if [ $singlescale ]; then
    rargs="$rargs --singlescale"
fi
if [ $singlecol ]; then
    rargs="$rargs --singlecol $singlecol"
fi
if [ $nogradient ]; then
    rargs="$rargs --nogradient"
fi
Rscript $scriptdir/blast_plot.R $rargs


