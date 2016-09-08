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

#3) Combine old and new reads
#mkdir combined_reads
#mkdir combined_reads/fastq

COMB_READ_DIR=combined_reads/fastq
OLD_READ_DIR=../HJ_D.tertiolecta/trimmed_reads/fastq
NEW_READ_DIR=trimmed_reads/fastq

#for sample in 'wt S1' 'L1 S2' 'L2 S3'
#do
 # set -- $sample
#  cat $OLD_READ_DIR/${1}_${2}_R1_trimmed.fastq.gz $NEW_READ_DIR/${1}_${2}_R1_trimmed.fastq.gz > $COMB_READ_DIR/${1}_${2}_R1_trimmed.fastq.gz
#  cat $OLD_READ_DIR/${1}_${2}_R2_trimmed.fastq.gz $NEW_READ_DIR/${1}_${2}_R2_trimmed.fastq.gz > $COMB_READ_DIR/${1}_${2}_R2_trimmed.fastq.gz
#done

#4) Double check genome size by doing kmer counting with Jellyfish
READ_DIR=$COMB_READ_DIR
JF_DIR=genome_size_estimate

#mkdir $JF_DIR

#for k in 17 21 25 31
#do
#  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/wt_kmer_counts_${k}.jf -C <(zcat $READ_DIR/wt_S1_R1_trimmed.fastq.gz)
#  jellyfish histo -o $JF_DIR/wt_kmer_counts_${k}.hist $JF_DIR/wt_kmer_counts_${k}.jf

#  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/L1_kmer_counts_${k}.jf -C <(zcat $READ_DIR/L1_S2_R1_trimmed.fastq.gz)
#  jellyfish histo -o $JF_DIR/L1_kmer_counts_${k}.hist $JF_DIR/L1_kmer_counts_${k}.jf

#  jellyfish count -m $k -s 150M --bf-size 3G -t 8 -o $JF_DIR/L2_kmer_counts_${k}.jf -C <(zcat $READ_DIR/L2_S3_R1_trimmed.fastq.gz)
#  jellyfish histo -o $JF_DIR/L2_kmer_counts_${k}.hist $JF_DIR/L2_kmer_counts_${k}.jf
#done

#Rscript $JF_DIR/kmer_hist_plot.R > $JF_DIR/JF_gensize_estimates.txt

#5) Another way of looking at coverage and estimating genome size - khmer read coverage plots: uses median kmer abund per read, so avoids issues with heterozygosity

source ~/Downloads/virtual_env/khmerEnv/bin/activate

#mkdir $JF_DIR/khmer

#for sample in 'wt S1' 'L1 S2' 'L2 S3'
#do
  #set -- $sample
  #load-into-counting.py -N 4 -x 8e8 -k 20 -T 8 $JF_DIR/khmer/${1}_R1.kh $READ_DIR/${1}_${2}_R1_trimmed.fastq.gz
  #count-median.py $JF_DIR/khmer/${1}_R1.kh $READ_DIR/${1}_${2}_R1_trimmed.fastq.gz $JF_DIR/khmer/${1}_R1_cov_per_read.txt
  #Rscript $JF_DIR/khmer/wt_cov_plot.R $JF_DIR/khmer/${1}_R1_cov_per_read.txt
#done > $JF_DIR/khmer/khmer_gensize_estimates.txt

#6) Check level of PhiX contamination

#mkdir phiX_screen

#fastq_screen --outdir phiX_screen ../combined_reads/fastq/*.fastq.gz
#grep PhiX phiX_screen/*.txt > phiX_screen/phiX_summary.txt

#Result: roughly 0.3% phiX contamination in each sample. Don't think it's worth the hassle to remove.

#7) WT assembly with SPAdes

READ_DIR=$COMB_READ_DIR

#mkdir assembly
#mkdir assembly/spades
#ASSEMBLY_DIR=assembly/spades

#spades.py -1 $READ_DIR/wt_S1_R1_trimmed.fastq.gz -2 $READ_DIR/wt_S1_R2_trimmed.fastq.gz -o $ASSEMBLY_DIR/Dt_wt_assembly --careful -t 8 > $ASSEMBLY_DIR/wt_spades_log.txt

#python ~/Downloads/quast-4.0/quast.py -o $ASSEMBLY_DIR/Dt_wt_assembly/quast -t 8 -s $ASSEMBLY_DIR/Dt_wt_assembly/scaffolds.fasta

#8) Assembly with Megahit
#mkdir assembly/megahit
#ASSEMBLY_DIR=assembly/megahit

#megahit -1 $READ_DIR/wt_S1_R1_trimmed.fastq.gz -2 $READ_DIR/wt_S1_R2_trimmed.fastq.gz --presets bulk -o $ASSEMBLY_DIR/Dt_wt_assembly > $ASSEMBLY_DIR/wt_spades_log.txt

#python ~/Downloads/quast-4.0/quast.py -o $ASSEMBLY_DIR/Dt_wt_assembly/quast -t 8 -s $ASSEMBLY_DIR/Dt_wt_assembly/final.contigs.fa


#9) Assembly with Megahit - reduce min-count

#mkdir assembly/megahit_mincount2
#ASSEMBLY_DIR=assembly/megahit_mincount2

#megahit -1 $READ_DIR/wt_S1_R1_trimmed.fastq.gz -2 $READ_DIR/wt_S1_R2_trimmed.fastq.gz --min-count 2 --k-min 21 --k-max 121 --k-step 10 --prune-level 2 -o $ASSEMBLY_DIR/Dt_wt_assembly

#10) Assembly with dipSPAdes

mkdir assembly/dipspades
ASSEMBLY_DIR=assembly/dipspades
dipspades.py -1 $READ_DIR/wt_S1_R1_trimmed.fastq.gz -2 $READ_DIR/wt_S1_R2_trimmed.fastq.gz -o $ASSEMBLY_DIR/Dt_wt_assembly --careful -t 8 > $ASSEMBLY_DIR/wt_dipspades_log.txt


#9) Compare assemblies

SPADES=assembly/spades/Dt_wt_assembly/scaffolds.fasta
MEGAHIT=assembly/megahit/Dt_wt_assembly/final.contigs.fa
MEGAHIT_MINCOUNT2=assembly/megahit_mincount2/Dt_wt_assembly/final.contigs.fa

#python ~/Downloads/quast-4.0/quast.py -o assembly/quast_comparison -t 8 -l SPAdes,Megahit,Megahit-mincount2 -s $SPADES $MEGAHIT $MEGAHIT_MINCOUNT2 