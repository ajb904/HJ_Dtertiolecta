#Start point: all relevant raw fastq files in raw_sequence/fastq directory.

#1) QC check of data using fastqc
#python scripts/QC_dir.py raw_sequence/fastq raw_sequence/QC

#2) Trim adapters from sequences, and quality trim
## NOTE: I upgraded to cutadapt 1.9.1 - previously 1.8.1
#mkdir trimmed_reads
#mkdir trimmed_reads/fastq
#mkdir trimmed_reads/QC
#python scripts/trim_fastq.py raw_sequence/fastq trimmed_reads/fastq TruSeq 20 50 > trim_log
#python scripts/QC_dir.py trimmed_reads/fastq trimmed_reads/QC

#2) Double check genome size by doing kmer counting with Jellyfish
READ_DIR=trimmed_reads/fastq
JF_DIR=genome_size_estimate

mkdir $JF_DIR

for k in 17 21 25 31
do
  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/wt_kmer_counts_${k}.jf -C <(zcat $READ_DIR/wt_S1_R1_trimmed.fastq.gz)
  jellyfish histo -o $JF_DIR/wt_kmer_counts_${k}.hist $JF_DIR/wt_kmer_counts_${k}.jf

  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/L1_kmer_counts_${k}.jf -C <(zcat $READ_DIR/L1_S2_R1_trimmed.fastq.gz)
  jellyfish histo -o $JF_DIR/L1_kmer_counts_${k}.hist $JF_DIR/L1_kmer_counts_${k}.jf

  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/L2_kmer_counts_${k}.jf -C <(zcat $READ_DIR/L2_S3_R1_trimmed.fastq.gz)
  jellyfish histo -o $JF_DIR/L2_kmer_counts_${k}.hist $JF_DIR/L2_kmer_counts_${k}.jf
done

Rscript $JF_DIR/kmer_hist_plot.R > $JF_DIR/JF_gensize_estimates.txt

#12) Another way of looking at coverage and estimating genome size - khmer read coverage plots: uses median kmer abund per read, so avoids issues with heterozygosity

source ~/Downloads/virtual_env/khmerEnv/bin/activate

mkdir $JF_DIR/khmer

for sample in 'wt S1' 'L1 S2' 'L2 S3'
do
  set -- $sample
  load-into-counting.py -N 4 -x 8e8 -k 20 -T 8 $JF_DIR/khmer/${1}_R1.kh $READ_DIR/${1}_${2}_R1_trimmed.fastq.gz
  count-median.py $JF_DIR/khmer/${1}_R1.kh $READ_DIR/${1}_${2}_R1_trimmed.fastq.gz $JF_DIR/khmer/${1}_R1_cov_per_read.txt
  Rscript $JF_DIR/khmer/wt_cov_plot.R $JF_DIR/khmer/${1}_R1_cov_per_read.txt
done > $JF_DIR/khmer/khmer_gensize_estimates.txt