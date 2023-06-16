#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:5                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=64G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/9.2.0
module load bedtools/2.30.0  

# prior step: split possorted_genome_bam.bam into sub bams to enable bedtools coverage without crashing due to memory limit
# 
# samtools view -bh possorted_genome_bam.bam chrM > chrM.bam
# samtools view -bh possorted_genome_bam.bam chr12 > chr12.bam
# samtools view -bh possorted_genome_bam.bam chr21 > chr21.bam

# calculate pseudobulk coverage for target regions from cellranger genome-aligned bam file output using sub bams
bedtools coverage -d -a ../targets.bed -b chrM.bam  > chrM.coverage.bed
bedtools coverage -d -a ../targets.bed -b chr12.bam  > chr12.coverage.bed
bedtools coverage -d -a ../targets.bed -b chr21.bam  > chr21.coverage.bed

# merge coverage
cat chrM.coverage.bed | grep chrM > coverage.bed
cat chr12.coverage.bed | grep chr12 >> coverage.bed
cat chr21.coverage.bed | grep chr21 >> coverage.bed
