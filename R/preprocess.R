getBins <- function(genome, chr=NULL, n=NULL, bin.size=NULL, offset=0) {
  stopifnot(!all(c(is.null(n), is.null(bin.size)), "specify either bin size or number of bins"))
  if (is.null(chr)) {
    chr = names(seqlengths(genome))
  }
  if (!is.null(n)) {
    bin.size = floor((seqlengths(genome)[chr] - offset) / n)
    names(bin.size) = chr
    n = rep(n, length(chr))
    names(n) = chr
  } else {
    n = floor((seqlengths(genome)[chr] - offset) / bin.size)
    names(n) = chr
    bin.size = rep(bin.size, length(chr))
    names(bin.size) = chr
  }
  
  g = GRanges()
  for (ch in chr) {
    g = c(g, GRanges(seqnames=ch, IRanges(start=0:(n[ch] - 1) * bin.size[ch] + 1 + offset, width=bin.size)))
  }
  return(g)
}

countBamInGRanges <- function(bam.file, granges, min.mapq=NULL, read.width=1) {
  require(GenomicRanges)
  require(Rsamtools)
  
  rds.counts <- numeric(length(granges))
  seq.names <- unique(as.character(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges)==seq.name]
      strand(granges.subset) <- "*"
      rds <- scanBam(bam.file,param=ScanBamParam(what=c("pos","mapq"),which=range(granges.subset)))
      if (!is.null(min.mapq)) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test = rep(T, length(rds[[1]]$mapq))
      }
      mapq.test = mapq.test & !is.na(rds[[1]]$pos)
      if (sum(mapq.test) > 0) {
        rds.ranges <- GRanges(seq.name,IRanges(start=rds[[1]]$pos[mapq.test],width=read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset,rds.ranges)
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- rds.counts.seq.name
      } else {
        rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
      }
    } else {
      rds.counts[as.logical(seqnames(granges)==seq.name)] <- 0
    }
  }
  rds.counts
}

countBedInGRanges <- function(bed.file, granges, read.width=1) {
  require(GenomicRanges)
  require(Rsamtools)

  # this is very memory inefficient but who cares!
  bed = bed2GR(bed.file)
  bed = resize(bed, read.width)
  
  return(countOverlaps(granges, bed))
}

extract <- function (pattern, string, perl = TRUE) {
  r <- paste(".*", pattern, ".*", sep = "")
  matched <- grep(r, string, perl = perl)
  result <- rep(NA, length(string))
  result[matched] <- sub(r, "\\1", string[matched], perl = perl)
  return(result)
}

gff2GR <- function(filename, gffAttrNames=NULL) {
  # read gff into genomic ranges
  require(GenomicRanges)
  regions = scan(filename, what=list(character(), character(), character(), numeric(), numeric(), character(), character(), character(), character()), comment.char="#", sep="\t")

  strand = regions[[7]]
  strand[strand == "."] = "*"
  strand[strand == "1"] = "+"
  strand[strand == "-1"] = "-"

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[4]], end=regions[[5]]),
    strand=strand)

  if (!is.null(gffAttrNames)) {
    df = data.frame(sapply(gffAttrNames, function(n)
      extract(paste(n , "=(.+?)(;|$)", sep=""), regions[[9]])),
      stringsAsFactors=F)
    elementMetadata(gr) = df
  }
  return(gr)
}

bed2GR <- function(filename, nfields=NULL) {
  # read bed into genomic ranges
  require(GenomicRanges)
  what = list(character(), numeric(), numeric(), character(), numeric(), character())
  if (is.null(nfields)) {
    hdr = read.csv(filename, sep="\t", header=F, nrows=3)
    nfields = ncol(hdr)
  }
  if (nfields > length(what)) {
    for (i in (length(what) + 1):nfields) {
      what[[i]] = character()
    }
  }
  regions = scan(filename, what=what[1:nfields])

  if (nfields >= 6) {
    strand = regions[[6]]
    strand[strand == "."] = "*"
  } else {
    strand = "*"
  }

  gr = GRanges(seqnames=regions[[1]],
    ranges=IRanges(start=regions[[2]], end=regions[[3]]),
    strand=strand)

  if (nfields >= 4) {
    names(gr) = regions[[4]]
    if (nfields >= 5) {
      m = DataFrame(score=as.numeric(regions[[5]]))
      if (nfields >= 7) {
        for (i in 7:nfields) {
          m = DataFrame(m, regions[[i]])
          colnames(m)[ncol(m)] = paste("field", i, sep="")
        }
      }
      elementMetadata(gr) = m
    }
  }
  
  return(gr)
}

preprocess.for.hmm <- function(signal.bam, genes, bin.size, genome, chr=NULL, expr=NULL, offset=0, n.expr.bins=5, plot=FALSE) {
  require(GenomicRanges)
  require(MASS)

  # define the binning
  bins = getBins(genome, chr=chr, bin.size=bin.size, offset=offset)
  
  if (!is.null(expr)) {
    ## get the subset of genes that have expression data
    expr.genes = genes[elementMetadata(genes)[,"ID"] %in% names(expr)]
    expr = expr[elementMetadata(expr.genes)[,"ID"]]
    ## create a binning for the expression values
    expr.bin.breaks = quantile(expr, prob=(0:(n.expr.bins - 1)) / n.expr.bins)
    expr.bin.fun = stepfun(expr.bin.breaks, 0:n.expr.bins)
    expr.bin = expr.bin.fun(expr[elementMetadata(expr.genes)[,"ID"]])
  } else if (!is.null(genes)) {
    expr.genes = genes
    expr.bin = rep(NA, length(genes))
  }
  
  # count the file
  if (length(grep(".bed$", signal.bam)) > 0) {
    signal = countBedInGRanges(signal.bam, bins)
  } else {
    signal = countBamInGRanges(signal.bam, bins)
  }

  bin.type = rep(NA, length(bins))
  bin.expr = rep(NA, length(bins))

  if (!is.null(genes)) {
    ov = findOverlaps(bins, expr.genes, type="within") # somehow as.matrix
    bins2genes = cbind(queryHits(ov), subjectHits(ov)) # does not work properly
    
    ## now we classify bins whether they overlap genes
    ## and if they do which range of expression they belong to
    bin.type = rep("intergenic", length(bins))
    bin.expr = rep(NA, length(bins))
    if (is.null(expr)) {
      bin.type[bins2genes[,1]] = "gene"
    } else {
      bin.type[bins2genes[,1]] = paste("expressed", expr.bin[bins2genes[,2]], sep=".")
      bin.expr[bins2genes[,1]] = expr[elementMetadata(expr.genes)[bins2genes[,2], "ID"]]
    }

    ## also flag bins that overlap genes
    ig = bin.type == "intergenic"
    bin.type[ig][queryHits(findOverlaps(bins[ig], expr.genes))] = "partial"
    bin.type = factor(bin.type)
  }
  
  tab = data.frame(chrom=as.character(seqnames(bins)),
      start=as.numeric(start(bins)), end=as.numeric(end(bins)),
      signal, bin.type, bin.expr)
  
  if (plot && !is.null(expr)) {
    par(mfrow=c(2, 1))
    hist(expr, breaks=50, main="Expression", xlab="Expression")
    abline(v=expr.bin.breaks, lty="dotted")
    boxplot(log10(signal + 1) ~ factor(bin.type, levels=c(paste("expressed", 1:5, sep="."), "intergenic", "partial"), labels=c(paste("Expr", 1:5, sep="."), "intergenic", "partial")), las=2, ylab="ChIP-seq", main="ChIP-seq")
  }
  
  return(tab)
}


makeGffAttributes <- function(df, cols=NULL) {
  if (ncol(df) == 0)
    return(rep("", nrow(df)))
  if (is.null(cols))
    cols = colnames(df)
  nvpairs = sapply(cols, function(s) paste(gsub(".", "_", s, fixed=T), df[,s], sep="="))
  if (!is.matrix(nvpairs)) {
    nvpairs = matrix(nvpairs, ncol=length(cols))
  }
  return(apply(nvpairs, 1, paste, collapse=";"))
}


GR2gff <- function(regions, filename, feature.type="experimental_feature", src="GenomicRanges", score=".", phase=".") {
  require(GenomicRanges)

  strnd = as.character(strand(regions))
  strnd[strnd == "*"] = "."
  
  tab = data.frame(as.character(seqnames(regions)), src, feature.type, as.numeric(start(regions)), as.numeric(end(regions)), score, strnd, phase, makeGffAttributes(as.data.frame(elementMetadata(regions))), stringsAsFactors=F)
  options(scipen=100)
  write.table(tab, file=filename, sep="\t", quote=F, row.names=F, col.names=F)
  options(scipen=0)
}


