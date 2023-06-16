#!/bin/bash

# example code to run souporcell on aligned bam file to extract donor and recipient using singularity container
# https://github.com/wheaton5/souporcell

#SBATCH -c 12                              # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                        # Partition to run in
#SBATCH --mem=16G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.souporcell.out      # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.souporcell.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL

GRCH38_REF=/home/lp175/GRCH38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
SNP_REF=/home/lp175/snp_reference/common_variants_grch38_with_chr.vcf
SOUPORCELL_CONTAINER=/n/app/singularity/containers/lp175/souporcell.sif  

singularity exec $SOUPORCELL_CONTAINER souporcell_pipeline.py -i possorted_genome_bam.bam -b barcodes.tsv -f $GRCH38_REF -t 12 -o souporcell -k 2 --common_variants $SNP_REF --skip_remap SKIP_REMAP
