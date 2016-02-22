#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description = 'Parse blast output to generate (binned) values for plotting, and/or summary values (reference positions covered in blast output and presence/absence)')
parser.add_argument('-i', '--infile', help = 'Input file (blast output in tab separated format)', required = True)
parser.add_argument('-c', '--contigfile', help = 'File containing list of contig names and sizes', required = True)
parser.add_argument('-o', '--binoutfile', help = 'File to write binned values to. If not specified, this step not performed')
parser.add_argument('-b', '--binsize', help = 'Window size for binning. Only relevant if --binoutfile is specified. Default is no binning.', type = int, default = 1)
parser.add_argument('-s', '--summ', help = 'Print summary values to stdout: number of sites present, percentage present, overall presence/absence, according to thresholds below', action = 'store_true')
parser.add_argument('--idthresh', help = 'Identity cutoff, default: 99 percent. Only relevant if --summ is set. Must be >0.', type = float, default = 99)
parser.add_argument('--lenthresh', help = 'Length (proportion of reference) cutoff, default: 80 percent. Only relevant if --summ is set', type = float, default = 80)
args = parser.parse_args()

contigs = {}  # key: contig name, value: list of the same length as this contig, containing %id of blast hit (if any) for each position in the contig (if no blast hit covering position of interest, 0; if multiple overlapping blast hits covering position of interest, maximum value)
if args.binoutfile: bins = {}  # Same as above, with values averaged over bin length (so each list has contig_length/bin_size elements). Note that any remaining contig length after last full bin is ignored - usually this simplification is inconsequential, however if there are many contigs and a large bin size is used, then it may cause artefacts where the reference length is apparently shortened.
clist = []  # list of contig names
totallength = 0  # combined length of all contigs
# Read in contig names and lengths to initialise the above
for line in open(args.contigfile):
    contig, length = line.split()
    length = int(length)
    totallength += length
    contigs[contig] = [0 for i in range(0, length)]  # Each element of each list initialised to 0 (default value to use if there are no blast hits covering this position)
    if args.binoutfile: bins[contig] = [0 for i in range(0, length / args.binsize)]  # Initialised to 0 as above
    clist.append(contig)

# For each blast hit, update the identity value for each position within the hit
for line in open(args.infile):
    hit = line.split()
    qcontig = hit[0]
    hitid = float(hit[2])
    hitlen = int(hit[3])
    qstart = int(hit[6])
    qend = int(hit[7])
    for pos in range(qstart, qend + 1):
        contigs[qcontig][pos-1] = max(contigs[qcontig][pos-1], hitid)

# If generating binned file, calculate average values within each bin and write to output file
if args.binoutfile:
    concatbins = []
    for c in clist:
        for i in range(0, len(contigs[c]) - (len(contigs[c]) % args.binsize)):
            bins[c][i / args.binsize] += contigs[c][i] / args.binsize
        concatbins += bins[c]
    with open(args.binoutfile, 'w') as binouthandle:
        binouthandle.write(','.join(map(str, concatbins)) + '\n')

# Generate summary
if args.summ:
    numsitespresent = len([idval for idlist in contigs.values() for idval in idlist if idval >= args.idthresh])
    proppresent = 100 * float(numsitespresent) / totallength
    overallpresence = 'Present' if proppresent >= args.lenthresh else 'Absent'
    print('\t'.join(map(str,[numsitespresent, '{0:.2f}'.format(proppresent), overallpresence])))


