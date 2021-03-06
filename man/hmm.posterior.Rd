\name{hmm.posterior}
\alias{hmm.posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  HMM posterior
}
\description{
  compute the posterior probability of a lognormal emission hmm
}
\usage{
hmm.posterior(x, meanlog, sdlog, max.it=500, eps=1e-4, rho=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    observation vector, can be also a two row matrix with bivariate
    observations
  }
  \item{meanlog}{
    vector holding the means of the lognormal emisson distributions, can
    be also a two row matrix with mean vectors in the columns
  }
  \item{sdlog}{
    vector holding the sd of the lognormal emisson distributions, can be
    also a two row matrix with the sqrt of the diagonal elements of sigma
  }
  \item{max.it}{
    maximum number of iterations
  }
  \item{eps}{
    stop iterations if the gain in log likelihood is less than \code{eps}
  }
  \item{rho}{
    if the bivariate model is used this is a vector of correlation
    coefficients for each state
  }
}
\details{
  posterior decoding
}
\value{
 matrix (T x N) of posterior probabilities
}
\references{
  Rabiner et al
}
\author{
  Matthias Heinig
}
\note{
}

\seealso{
\code{\link{zinba.hmm.posterior}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, logmean, logsd, A, p_i) {return(posterior)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
