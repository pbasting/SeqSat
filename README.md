# SeqSat
produces a saturation curve for RNAseq datasets

### Dependencies
* BWA
* SAMtools
* BEDtools

### Installation
* download anaconda from : https://www.anaconda.com/distribution/
* create environment using the provided `environment.yml` file
```
conda env create -f environment.yml
conda activate seqSat
```
* download SeqSat repository
```
git clone git@github.com:pbasting/SeqSat.git
```
### Usage
```
python3 seqSat.py \
    -i SRR111111.fq.gz SRR2222222.fq.gz \
    -f reference.fasta \
    -g genes.gff \
    -n sample1 sample2 \
    -o plot_dir
```
