\name{run.univariate.hmm}
\alias{run.univariate.hmm}
\title{
Univariate HMM analysis
}
\description{
This function provides a highlevel interface to the univariate HMM for
single sample analysis of broad histone marks. It takes the count data
as input, estimates the model parameters and outputs the posterior
probabilities of the model.
}
\usage{
run.univariate.hmm(fname, data = NULL, n.expr.bins = 5, em = FALSE, chrom = NULL, maxq = 1 - 1e-04, redo = FALSE, model.constructor = "Zinba", baum.welch = FALSE)
}

\arguments{
  \item{fname}{
    filename of the count data. This should end with ".txt" since it is
    also used to derive the output file names. The file should have the
    format produced by \code{\link{preprocess.for.hmm}}.
}
  \item{data}{
    if the data is already in memory, it can be passed as an argument
    and no input file will be read. Note that you still have to provide
    the \code{fname} argument as this will be used to derive output filenames.
}
  \item{n.expr.bins}{
    when using expression data for the parameter estimation provide the
    number of expression bins (same as used in
    \code{\link{preprocess.for.hmm}}).  
}
  \item{em}{
    \code{logical} indicating whether to use the EM algorithm for
    parameter estimation.
}
  \item{chrom}{
    list of chromosomes to use for the parameter estimation.
}
  \item{maxq}{
    maximum quantile of counts. This is used to achieve better numerical
    stability and to remove mapping artifacts with very high
    counts. Specifically this works by determining the quantile \code{k
      = quantile(counts, maxq)} and then setting \code{counts[counts > k] = k}.
}
  \item{redo}{
    logical indicating whether to redo the parameter estimation and
    posterior computation. If FALSE the precomputed results will be loaded.
}
  \item{model.constructor}{
    character vector indicating the type of distributions to be used in
    the EM algorithm. This can be one of "Zinba", "Negbinom",
    "LogNormal", "Zero". If length is one, the same distribution will be
    used for all components.
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
  
  Some notes on using the EM algorithm. When you plan to fit a two
  component model with \code{model.constructor=c("Zinba", "Negbinom")}
  it is in practice better to fit a three component model with
  \code{model.constructor=c("Zero", "Negbinom", "Negbinom")}. This is
  equivalent since the zero inflated negative binomial component is
  itself a mixture of a zero component and a negative binomial.
}
\value{
  \code{data.frame} with the original input data and the posterior
  probabilities. This is returned invisibly.
}
\references{http://histonehmm.molgen.mpg.de}
\author{
Matthias Heinig
}
\note{
%%  ~~further notes~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(rat.H3K27me3)
posterior = run.univariate.hmm("BN_H3K27me3.txt", data=BN, n.expr.bins=5, em=TRUE, chrom="chr19", maxq=1-1e-04, redo=TRUE, model.constructor="Zinba", baum.welch=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
