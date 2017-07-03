#!/usr/bin/env Rscript

library(optparse)
library(histoneHMM)
library(GenomicRanges)

## parse command line options
option_list <- list(
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help=paste("gene annotation in gff/gtf format with ",
                  "ID/gene_id field (optinal)", collapse=""))
                    ,
  make_option(c("-b", "--binsize"), type="numeric", default=1000,
              help="binsize for counting reads [default \"%default\"]"),
  make_option(c("-c", "--chromlen"), type="character", default=NULL,
              help="file with chromosome lengths (<chr><tab><len>) no header"),
  make_option(c("-e", "--expression"), default=NULL,
              help = paste("gene expression values (optional) format:",
                  "(<gene_id><tab><expression_value>) no header", collapse="")),
  make_option("--em", type="logical", default=TRUE,
              help=paste("use the EM algorithm to fit densities, if FALSE",
                  "gene expression data is used and has to be provided",
                  "[default \"%default\"]", collapse="")),
  make_option(c("-o", "--outprefix"), default=NULL,
              help="prefix for all output files [default \"histoneHMM\"]"),
  make_option(c("-p", "--preprocessed"), default=NULL,
              help="file with previously preprocessed data"),
  make_option(c("-t", "--train"), default=NULL,
              help="chromosome used for training the model (optional)"),    
  make_option("--verbose", action="store_true", default=FALSE,
              help="more output"),
  make_option(c("-P", "--probability"), default=0.5,
              help="threhshold for the posterior probability to call regions")
)

parser = OptionParser(usage="usage: %prog [options] bamfile",
  option_list=option_list)
opt <- parse_args(parser, positional_arguments=T)
if (opt$options$verbose) {
  print(opt)
}

## check if we need to do the preprocessing
if (is.null(opt$options$preprocessed)) {

  ## check that we have all input data for the preprocessing
  if (length(opt$args) < 1) {
    print_help(parser)
    q("no")
  }
  
  bamfile = opt$args[1]
  prefix = opt$options$outprefix
  if (is.null(prefix)) {
    prefix = "histoneHMM"
  }
  
  if (is.null(opt$options$chromlen)) {
    cat("Please specify the chromlen file (mandatory)\n")
    print_help(parser)
    q("no")
  } else {
    if (!file.exists(opt$options$chromlen)) {
      cat(paste("chromlen file '", opt$options$chromlen,
                "' not found\n", sep=""))
      q("no")
    }
  }
  
  chroms = read.table(opt$options$chromlen, stringsAsFactors=FALSE)
  sl = chroms[,2]
  names(sl) = chroms[,1]
  genome = GRanges(seqlengths=sl)
  
  if (!is.null(opt$options$annotation)) {
    ## check if we have gff or gtf format
    gff = grep("gff$", opt$options$annotation)
    if (length(gff) > 0) {
      genes = gff2GR(opt$options$annotation, "ID")
    } else {
      genes = gtf2GR(opt$options$annotation, "gene_id")
    }
  } else {
    genes = NULL
  }
  
  if (!is.null(opt$options$expression)) {
    if (!file.exists(opt$options$expression)) {
      cat(paste("gene expression file '", opt$options$expression,
                "' not found\n", collapse=""))
      q("no")
    }
    if (is.null(genes)) {
      cat("gene annotation not specified! (mandatory if expression is given)\n")
      print_help(parser)
      q("no")
    }
    expr = read.table(opt$options$expression, row.names=1)
    n = rownames(expr)
    expr = expr[,1]
    names(expr) = n

    ## open a device for the boxplot generated while preprocessing
    ## (only when expression is given)
    svg(file=paste(prefix, "-expr-boxplot.svg", sep=""))
  } else {
    expr = NULL
  }
  data = preprocess.for.hmm(bamfile, genes, bin.size=opt$options$binsize,
      genome, expr=expr, plot=T)
  write.table(data, paste(prefix, ".txt", sep=""), sep="\t", quote=F,
              row.names=F)
  if (!is.null(expr)) {
    dev.off()
  }
} else {
  if (!file.exists(opt$options$preprocessed)) {
    cat(paste("preprocessed file '", opt$options$preprocessed,
              "' not found\n", sep=""))
    q("no")
  }
  ## when no prefix is given, we use the name of the preprocessed file
  if (is.null(opt$options$outprefix)) {
    prefix = strsplit(basename(opt$options$preprocessed), ".", fixed=T)[[1]]
    n = length(prefix) - 1
    if (n == 0) {
      n = 1
    }
    prefix = paste(prefix[1:n], collapse=".")
    prefix = file.path(dirname(opt$options$preprocessed), prefix)
  } else {
    prefix = opt$options$outprefix
  }
  data = read.csv(opt$options$preprocessed, sep="\t", stringsAsFactors=F)
}


## do the actual analysis
posterior = run.univariate.hmm(paste(prefix, ".txt", sep=""), data=data,
    em=opt$options$em, chrom=opt$options$train, maxq=1-1e-3,
    model.constructor=c("Zero", "Negbinom", "Negbinom"))

regions = callRegions(posterior, opt$options$probability,
  posterior.col="lowly.expressed", expr.col=NULL)

GR2gff(regions, paste(prefix, "-regions.gff", sep=""))
