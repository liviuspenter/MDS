#!/bin/bash

# code run on O2 cluster to generate pseudoaligments of scRNA-seq libraries to quantify PCR handles

python /Users/shaka87/dfci/github/kite/featuremap/featuremap.py ../../../ref/seqs.csv  --header

kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa

kallisto bus -i FeaturesMismatch.idx -o ./ -x 10xv3 -t 12 bamtofastq_S1_L001_R* bamtofastq_S1_L002_R* bamtofastq_S1_L003_R* bamtofastq_S1_L004_R* 

bustools sort -t 4 -o ./output_sorted.bus ./output.bus
mkdir ./featurecounts/
bustools count -o featurecounts/featurecounts --genecounts -g ./FeaturesMismatch.t2g -e ./matrix.ec -t ./transcripts.txt ./output_sorted.bus

rm output.bus
rm output_sorted.bus
