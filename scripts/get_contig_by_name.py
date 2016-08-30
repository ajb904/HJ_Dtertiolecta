#!/usr/local/bin/python

#Extract fasta contigs by name

import sys
from Bio import SeqIO

contig_wanted = sys.argv[1]

if sys.argv[2] == '-':
	fasta_file = sys.stdin
else:
	fasta_file = sys.argv[2]

#Read through fastq file and output reads with names that are in reads_wanted
for seq in SeqIO.parse(fasta_file, 'fasta'):
	if seq.id == contig_wanted:
		SeqIO.write(seq, sys.stdout, 'fasta')
