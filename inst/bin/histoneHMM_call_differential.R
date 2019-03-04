#!/usr/bin/env Rscript

library(optparse)
library(histoneHMM)
library(GenomicRanges)

## parse command line options
option_list <- list(
  make_option("--em", type="logical", default=TRUE,
              help=paste("use the EM algorithm to fit densities, if FALSE",
                  "gene expression data is used and has to be provided",
                  "[default \"%default\"]", collapse="")),
  make_option(c("-o", "--outdir"), default="histoneHMM_differential",
              help=paste("outputdir for all output files",
                "[default \"histoneHMM_differential\"]")),
  make_option("--sample1", default="sample1",
              help="name of sample1 [default \"sample1\"]"),
  make_option("--sample2", default="sample2",
              help="name of sample2 [default \"sample2\"]"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="more output"),
  make_option("--baumwelch", action="store_false",
              help="Whether the Baum Welch algorithm should be used for
                    parameter estimation.")
  make_option(c("-P", "--probability"), default=0.5,
              help="threhshold for the posterior probability to call regions")
)

parser = OptionParser(
  usage=paste("usage: %prog [options] preprocessed_sample1.txt",
    "preprocessed_sample2.txt"),
  option_list=option_list)
opt <- parse_args(parser, positional_arguments=T)
if (opt$options$verbose) {
  print(opt)
}

if (length(opt$args) < 2) {
  cat("Not enough arguments!\n")
  print_help(parser)
  q("no")
}

fname1 = opt$args[1]
fname2 = opt$args[2]

## run the HMM
regions = get.differential.regions(fname1, fname2, outdir=opt$options$outdir,
  sample1=opt$options$sample1, sample2=opt$options$sample2,
  em=opt$options$em, maxq=1-1e-3, baum.welch=opt$options$baumwelch,
  cutoff=opt$options$probability)
