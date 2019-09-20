# SeqSat
produces a saturation curve for RNAseq datasets

### Dependencies
* BWA
* SAMtools
* BEDtools

### Usage
```
python3 seqSat.py \
    -i SRR111111.fq.gz SRR2222222.fq.gz \
    -f reference.fasta \
    -g genes.gff \
    -n sample1 sample2 \
    -o plot_dir
```
