\name{em}
\alias{em}
\alias{EMResult-class}
\title{
Expectation maximization algorithm
}
\description{
Expectation maximization algorithm for different classes of
distributions. Currently only the Zinba distribution is used, but
Lognormal is also available.
}
\usage{
em(x, ncomp, prop = NULL, maxit = 100, eps = 1e-04, model.constructor = "Zinba")
}
\arguments{
  \item{x}{
data vector
}
  \item{ncomp}{
number of mixture components
}
  \item{prop}{
initial guess of the mixture proportions
}
  \item{maxit}{
maximum number of iterations
}
  \item{eps}{
stop if the gain in loglikelihood is less than \code{eps}
}
\item{model.constructor}{
name of a constructor for the mixture components, can be "Zinba" or "Lognormal"
}
}
\details{
  If \code{prop} is NULL a uniform distribution is assumed.
}
\value{
  returns an object of class EMResult which has the following slots:
  \item{models}{a list containing the models for each component}
  \item{logLikelihood}{a matrix with the logLikelihood for each data
    point (rows) and each mixture component (cols) }
  \item{dataLogLikelihood}{the total log likelihood of the data}
  \item{proportions}{the mixture proportions}
}
\references{
}
\author{
Matthias Heinig
}
\note{
}

\seealso{
The models are objects of the class \code{\linkS4class{Distribution}} for
example \code{\linkS4class{Zinba}} or \code{\linkS4class{Lognormal}}.
}
\examples{
x = rzinba(1000, size=10, mu=5, beta=0.1)
y = rzinba(1000, size=10, mu=10, beta=0.02)

emfit = em(c(x, y), ncomp=2, model.constructor="Zinba")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

