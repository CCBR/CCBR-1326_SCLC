# CCBR 1326

integrate ATAC-seq and RNA-seq data in small cell lung cancer

## Usage

on Biowulf

### integrate ATAC & RNAseq DGE

```sh
module load nextflow
nextflow run main.nf -resume
```

### tobias TF footprinting

```sh
bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/sovacool-dev-tobias/run_tobias.bash \
  -w=/data/SCLCgenomics/sovacoolkl/CCBR-1326/ccbr_tobias \
  -m=dryrun
```
