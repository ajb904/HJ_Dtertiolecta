#!/usr/local/bin/python

# Take paired end fastq files, and trim each using cutadapt
# Input: 1) input directory, 2) output directory, 3) TruSeq or Nextera, 4) quality threshold, 5) minimum length to keep

import sys
import os
import subprocess

in_dir = sys.argv[1]
out_dir = sys.argv[2]
lib_type = str(sys.argv[3]).lower()
qual = int(sys.argv[4])
min_len = int(sys.argv[5])

#Check library type, and use that to assign adapter sequences to trim
if lib_type == 'truseq':
    adapter_f = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter_r = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
elif lib_type == 'nextera':
    adapter_f = 'CTGTCTCTTATACACATCT'
    adapter_r = 'CTGTCTCTTATACACATCT'
else:
    raise AttributeError('Not a valid library type')


#Get names of all R1 fastq files from input directory
fq_basenames = [os.path.basename(f) for f in os.listdir(in_dir) if f.endswith('R1_001.fastq.gz')]

for f in fq_basenames:
    f1 = f
    f2 = f.replace('R1', 'R2')

    r1_in = os.path.join(in_dir, f1)
    r2_in = os.path.join(in_dir, f2)

    r1_out = os.path.join(out_dir, f1.replace('L001_R1_001','R1_trimmed'))
    r2_out = os.path.join(out_dir, f2.replace('L001_R2_001', 'R2_trimmed'))

    cline_args = ['-a', adapter_f,
                  '-A', adapter_r,
                  '-q', str(qual),
                  '-m', str(min_len),
                  '-o', r1_out,
                  '-p', r2_out]

    cline = ['cutadapt'] + cline_args + [r1_in, r2_in]

    subprocess.call(cline)
