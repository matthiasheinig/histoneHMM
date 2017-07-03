# histoneHMM

histoneHMM is a software to analyse ChIP-seq data of histone modifications with broad genomic footprints like H3K27me3. It allows for calling modified regions in single samples as well as for calling differentially modified regions in a comparison of two samples.

## News

Version 1.7 is an update that works with the new Rcpp package (Verion > 0.12)

Version 1.6 fixes some smaller bugs and removes the dependency on the GNU scientific library

We released version 1.5 introduces a command line interface (see vignette), improved preprocessing and an updated vignette.

We have just released the new and improved version 1.4 with a high level interface for more convenience.

## Installation

```{r}
install.packages("devtools")
devtools::install_github("matthiasheinig/histoneHMM")
```

## Usage

We provide a small test data set to try out the package:

### Preprocessed data:
[BN.txt](http://histonehmm.molgen.mpg.de/data/BN.txt)
[SHR.txt](http://histonehmm.molgen.mpg.de/data/SHR.txt)

### Raw data:
[BN.bam](http://histonehmm.molgen.mpg.de/data/BN.bam)
[BN.bam.bai](http://histonehmm.molgen.mpg.de/data/BN.bam.bai)
[SHR.bam](http://histonehmm.molgen.mpg.de/data/SHR.bam)
[SHR.bam.bai](http://histonehmm.molgen.mpg.de/data/SHR.bam.bai)
[chroms.txt](http://histonehmm.molgen.mpg.de/data/chroms.txt) : table of chromosome lengths
[ensembl59-genes.gff](http://histonehmm.molgen.mpg.de/data/ensembl59-genes.gff) : gene annotation in gff format (only gene bodies, no exon information)
[expression.txt](http://histonehmm.molgen.mpg.de/data/expression.txt): RNA-seq gene expression

## Highlevel interface

We have a new and improved version with a high level interface for more convenience. Please check out the package vignette:

[histoneHMM.pdf](http://histonehmm.molgen.mpg.de/v1.6/histoneHMM.pdf)
[histoneHMM.R](http://histonehmm.molgen.mpg.de/v1.6/histoneHMM.R)

## Citation

Heinig M, Colome-Tatche M, Rintisch C, Schäfer S, Pravenec M, Hubner N, Vingron M, Johannes F. histoneHMM: Differential analysis of histone modifications with broad genomic footprints. BMC Bioinformatics 2015; 16:60

## Contact

Matthias Heinig: matthias.heinig@helmholtz-muenchen.de