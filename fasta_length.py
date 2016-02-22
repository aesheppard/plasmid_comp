#!/usr/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = 'Returns the length of sequence in a fasta file.')
parser.add_argument('-f', '--fasta', help = 'Fasta file', required = True)
parser.add_argument('-t', '--total', help = 'Print total length only, instead of each contig individually', action = 'store_true')
parser.add_argument('-n', '--names', help = 'Print contig names as well as lengths. Ignored if --total is True.', action = 'store_true')
args = parser.parse_args()

if args.total:
    length = 0
    for contig in SeqIO.parse(args.fasta, 'fasta'):
        length += len(contig)
    print(length)

elif args.names:
    for contig in SeqIO.parse(args.fasta, 'fasta'):
        print('\t'.join([contig.id, str(len(contig))]))

else:
    for contig in SeqIO.parse(args.fasta, 'fasta'):
        print(len(contig))
