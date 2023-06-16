This directory contains the pseudoaligned counts generated using [kallisto | bustools](https://github.com/pachterlab/kallistobustools).

The fastq files of the scRNA-seq libraries were processed following standard processing steps:

### 1) Generation of FeatureMap using [reference](ref/seqs.csv).

```sh
python featuremap.py ./ref/seqs.csv  --header
```

### 2) Running kallisto

```sh
kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa

kallisto bus -i FeaturesMismatch.idx -o ./ -x 10xv3 -t 12 <fastq files>
```

### 3) Running bustools

```sh
bustools sort -t 4 -o ./output_sorted.bus ./output.bus
mkdir ./featurecounts/
bustools count -o featurecounts/featurecounts --genecounts -g ./FeaturesMismatch.t2g -e ./matrix.ec -t ./transcripts.txt ./output_sorted.bus
```
