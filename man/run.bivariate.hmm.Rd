\name{run.bivariate.hmm}
\alias{run.bivariate.hmm}
\title{Bivariate HMM analysis}
\description{
This function provides a highlevel interface to the bivariate HMM for
differential analysis of broad histone marks. It takes two count data sets
as input, estimates the model parameters and outputs the posterior
probabilities of the model.
}
\usage{
run.bivariate.hmm(fname1, fname2, outdir, data1 = NULL, data2 = NULL, sample1 = NULL, sample2 = NULL, n.expr.bins = 5, maxq = 1 - 1e-04, em = FALSE, chrom = NULL, baum.welch = FALSE)
}

\arguments{
  \item{fname1}{
    filename of the count data for sample1. The file should have the
    format produced by \code{\link{preprocess.for.hmm}}.
}
  \item{fname2}{
    filename of the count data for sample1. The file should have the
    format produced by \code{\link{preprocess.for.hmm}}.
}
  \item{outdir}{
    path to the output directory.
}
  \item{data1}{
    count data set for sample1. If the data is already in memory, it can
    be passed as an argument and no input file will be read.
}
  \item{data2}{
    count data set for sample2.
}
  \item{sample1}{
    name of sample1.
}
  \item{sample2}{
    name of sample2.
}
  \item{n.expr.bins}{
    when using expression data for the parameter estimation provide the
    number of expression bins (same as used in
    \code{\link{preprocess.for.hmm}}). 
}
  \item{maxq}{
    maximum quantile of counts. This is used to achieve better numerical
    stability and to remove mapping artifacts with very high
    counts. Specifically this works by determining the quantile \code{k
      = quantile(counts, maxq)} and then setting \code{counts[counts > k] = k}.
}
  \item{em}{
    \code{logical} indicating whether to use the EM algorithm for
    parameter estimation.
}
  \item{chrom}{
    list of chromosomes to use for the parameter estimation.
}
  \item{baum.welch}{
    logical indicating if the Baum Welch algorithm should be used for
    parameter estimation.
}
}
\details{
  As a side effect this function produces a number of files for the
  model parameters, plots that show the parameter fit and the posterior
  probabilities. 
}
\value{
  \code{data.frame} with the original input data and the posterior
  probabilities. This is returned invisibly.
}
\references{http://histonehmm.molgen.mpg.de
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fname1, fname2, outdir, data1 = NULL, data2 = NULL, 
    sample1 = NULL, sample2 = NULL, n.expr.bins = 5, maxq = 1 - 
        1e-04, em = FALSE, chrom = NULL, baum.welch = FALSE) 
{
    if (is.null(data1)) {
        data1 = read.csv(fname1, sep = "\t", stringsAsFactors = F)
    }
    if (is.null(data2)) {
        data2 = read.csv(fname2, sep = "\t", stringsAsFactors = F)
    }
    stopifnot(nrow(data1) == nrow(data2))
    max1 = round(quantile(data1$signal, p = maxq))
    data1$signal[data1$signal > max1] = max1
    max2 = round(quantile(data2$signal, p = maxq))
    data2$signal[data2$signal > max2] = max2
    if (is.null(chrom)) {
        chrom = unique(data1$chrom)
    }
    legacy = FALSE
    if (legacy) {
        both.high = which(data1$bin.type == paste("expressed", 
            n.expr.bins, sep = ".") & data2$bin.type == paste("expressed", 
            n.expr.bins, sep = "."))
        both.low = which(data1$bin.type == "expressed.1" & data2$bin.type == 
            "expressed.1")
        not.expressed = paste("expressed", 1:2, sep = ".")
        really.expressed = paste("expressed", 3:5, sep = ".")
        up = which(data1$bin.type \%in\% not.expressed & data2$bin.type \%in\% 
            really.expressed)
        down = which(data1$bin.type \%in\% really.expressed & data2$bin.type \%in\% 
            not.expressed)
    }
    else {
        posterior1 = run.univariate.hmm(fname1, data = data1, 
            n.expr.bins = n.expr.bins, em = em, chrom = chrom, 
            maxq = maxq, redo = FALSE, baum.welch = baum.welch)
        posterior2 = run.univariate.hmm(fname2, data = data2, 
            n.expr.bins = n.expr.bins, em = em, chrom = chrom, 
            maxq = maxq, redo = FALSE, baum.welch = baum.welch)
        threshold = 0.9
        while (threshold > 0.5) {
            cols = c("highly.expressed", "lowly.expressed")
            map1 = cols[apply(posterior1[, cols], 1, function(x) {
                gt = x > threshold
                if (any(gt, na.rm = T)) 
                  return(which(gt))
                else return(NA)
            })]
            map2 = cols[apply(posterior2[, cols], 1, function(x) {
                gt = x > threshold
                if (any(gt, na.rm = T)) 
                  return(which(gt))
                else return(NA)
            })]
            both.high = which(map1 == "highly.expressed" & map2 == 
                "highly.expressed")
            both.low = which(map1 == "lowly.expressed" & map2 == 
                "lowly.expressed")
            up = which(map1 == "lowly.expressed" & map2 == "highly.expressed")
            down = which(map1 == "highly.expressed" & map2 == 
                "lowly.expressed")
            if (length(up) < 50 || length(down) < 50) {
                threshold = threshold - 0.1
            }
            else {
                break
            }
        }
    }
    fit.low = fit.zinba.copula(data1$signal[both.low], data2$signal[both.low])
    fit.high = fit.zinba.copula(data1$signal[both.high], data2$signal[both.high])
    fit.lowx = fit.zinba.copula(data1$signal[up], data2$signal[up], 
        marginal.x = fit.low$marginal.x, marginal.y = fit.high$marginal.y)
    fit.lowy = fit.zinba.copula(data1$signal[down], data2$signal[down], 
        marginal.x = fit.high$marginal.x, marginal.y = fit.low$marginal.y)
    fits = list(unmod.both = fit.high, mod.x = fit.lowx, mod.y = fit.lowy, 
        mod.both = fit.low)
    zinba.params = sapply(fits, fit.to.col)
    dir.create(outdir)
    if (is.null(sample1)) {
        sample1 = "sample1"
    }
    if (is.null(sample2)) {
        sample2 = "sample2"
    }
    outname = paste(sample1, "-vs-", sample2, sep = "")
    write.table(zinba.params, file = file.path(outdir, paste(outname, 
        "-zinbacopula-params.txt", sep = "")), sep = "\t", quote = F)
    for (chrom in unique(data1$chrom)) {
        cat(chrom, "\n")
        this.chrom = data1$chrom == chrom
        signal1 = data1[this.chrom, "signal"]
        signal2 = data2[this.chrom, "signal"]
        signal = rbind(signal1, signal2)
        posterior.zinba = zinbacopula.hmm.posterior(signal, fits)
        colnames(posterior.zinba) = colnames(zinba.params)
        colnames(posterior.zinba) = gsub("mod.x", sample1, colnames(posterior.zinba))
        colnames(posterior.zinba) = gsub("mod.y", sample2, colnames(posterior.zinba))
        map.zinba = colnames(posterior.zinba)[apply(posterior.zinba, 
            1, which.max)]
        d1 = data1[this.chrom, c("signal", "bin.expr")]
        colnames(d1) = paste(colnames(d1), sample1, sep = ".")
        d2 = data2[this.chrom, c("signal", "bin.expr")]
        colnames(d2) = paste(colnames(d2), sample2, sep = ".")
        res = cbind(data1[this.chrom, c("chrom", "start", "end", 
            "bin.type")], d1, d2, posterior.zinba, map.zinba)
        first = (chrom == unique(data1$chrom)[1])
        write.table(res, file.path(outdir, paste(outname, ".txt", 
            sep = "")), sep = "\t", quote = F, row.names = F, 
            col.names = first, append = !first)
    }
    bivariate.posterior = read.csv(file.path(outdir, paste(outname, 
        ".txt", sep = "")), sep = "\t", stringsAsFactors = F)
    return(invisible(bivariate.posterior))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
