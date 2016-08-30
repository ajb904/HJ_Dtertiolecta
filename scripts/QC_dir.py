#!/usr/local/bin/python

import sys
import os
import subprocess

# Usage: QC_dir.py <fastq dir> <output dir> If output dir not given, use current dir

fastq_dir = sys.argv[1]
try:
    out_dir = sys.argv[2]
except IndexError:
    out_dir = os.curdir

cline = 'fastqc -o %s %s/*.fastq.gz' % (out_dir, fastq_dir)

subprocess.call(cline, shell=True)