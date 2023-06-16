#!/bin/bash

# code run on O2 cluster to generate pseudoaligments of scRNA-seq libraries to quantify PCR handles

#SBATCH -c 12                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                           # Partition to run in
#SBATCH --mem=16G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0
module load kallisto/0.46.2

#python /Users/shaka87/dfci/github/kite/featuremap/featuremap.py ../../../ref/seqs.csv  --header

#kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa

kallisto bus -i FeaturesMismatch.idx -o ./ -x 10xv3 -t 12 \
bamtofastq_S1_L001_R1_001.fastq.gz bamtofastq_S1_L001_R2_001.fastq.gz \
bamtofastq_S1_L001_R1_002.fastq.gz bamtofastq_S1_L001_R2_002.fastq.gz \
bamtofastq_S1_L001_R1_003.fastq.gz bamtofastq_S1_L001_R2_003.fastq.gz \
bamtofastq_S1_L002_R1_001.fastq.gz bamtofastq_S1_L002_R2_001.fastq.gz \
bamtofastq_S1_L002_R1_002.fastq.gz bamtofastq_S1_L002_R2_002.fastq.gz \
bamtofastq_S1_L002_R1_003.fastq.gz bamtofastq_S1_L002_R2_003.fastq.gz \
bamtofastq_S1_L003_R1_001.fastq.gz bamtofastq_S1_L003_R2_001.fastq.gz \
bamtofastq_S1_L003_R1_002.fastq.gz bamtofastq_S1_L003_R2_002.fastq.gz \
bamtofastq_S1_L003_R1_003.fastq.gz bamtofastq_S1_L003_R2_003.fastq.gz \
bamtofastq_S1_L004_R1_001.fastq.gz bamtofastq_S1_L004_R2_001.fastq.gz \
bamtofastq_S1_L004_R1_002.fastq.gz bamtofastq_S1_L004_R2_002.fastq.gz \
bamtofastq_S1_L004_R1_003.fastq.gz bamtofastq_S1_L004_R2_003.fastq.gz 

/home/lp175/software/bustools/bustools sort -t 4 -o ./output_sorted.bus ./output.bus
mkdir ./featurecounts/
/home/lp175/software/bustools/bustools count -o featurecounts/featurecounts --genecounts -g ./FeaturesMismatch.t2g -e ./matrix.ec -t ./transcripts.txt ./output_sorted.bus

rm output.bus
rm output_sorted.bus
