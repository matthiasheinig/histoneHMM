hmm.posterior <- function(x, meanlog, sdlog, max.it=500, eps=1e-4, rho=0) {

  # number of states
  N = length(meanlog)
  stopifnot(length(sdlog) == N)
  # length of observation sequence
  T = length(x)
  
  # we also allow bivariate input if x is a 2-row matrix
  bivariate = is.matrix(x) && nrow(x) == 2
  if (bivariate) {
    # make sure that meanlog and sdlog are also 2-row matrices
    stopifnot(is.matrix(meanlog), nrow(meanlog) == 2,
              is.matrix(sdlog), nrow(sdlog) == 2)
    N = ncol(meanlog)
    T = ncol(x)
    if (length(rho) == 1) {
      rho = rep(rho, N)
    }
    stopifnot(length(rho) == N)
  }
  
  # void R_hmm_posterior(double* O, int* T, int* N, double* mu, double* sigma, int* iterationMAX, double* eps, double* post)
  z = .C("R_hmm_posterior",
    as.double(x + 1), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(meanlog), # double* mu
    as.double(sdlog), # double* sigma
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    as.integer(bivariate), # int* bivariate
    double(length=N*N), # double* A
    double(length=N), # double* proba
    as.double(rho)) # double* rho
  # POST will hold the results

  posterior = matrix(z[[8]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[10]], ncol=N), proba=z[[11]])
  attr(posterior, "params") = params
  return(posterior)
}

zinba.hmm.posterior <- function(x, mu, size, beta, max.it=500, eps=1e-4) {

  # number of states
  N = length(mu)
  stopifnot(length(size) == N)
  stopifnot(length(beta) == N)
  # length of observation sequence
  T = length(x) 
 
  # void R_hmm_posterior_zinba(double* O, int* T, int* N, double* mu, double* size, double* beta, int* iterationMAX, double* eps, double* post, double* A, double* proba)
  z = .C("R_hmm_posterior_zinba",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(mu), # double* mu
    as.double(size), # double* size
    as.double(beta), # double* beta
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N)) # double* proba
  # POST will hold the results

  posterior = matrix(z[[9]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[10]], ncol=N), proba=z[[11]])
  attr(posterior, "params") = params
  return(posterior)
}


bivariate.zinba.hmm.posterior <- function(x, marginal.mu, marginal.size, marginal.beta, mu.coef, size.coef, beta0, max.it=500, eps=1e-4) {

  # number of states
  N = length(marginal.mu)
  stopifnot(length(marginal.size) == N)
  stopifnot(length(marginal.beta) == N)
  stopifnot(ncol(mu.coef) == N)
  stopifnot(ncol(size.coef) == N)
  stopifnot(length(beta0) == N)
  
  # length of observation sequence
  T = ncol(x) 
 
  #  void R_hmm_posterior_bivariatezinba(double* O, int* T, int* N, double* marginal_mu, double* marginal_size, double* marginal_beta, double* size_coef, double* mu_coef, double* beta0, int* iterationMAX, double* eps, double* post, double* A, double* proba)
  z = .C("R_hmm_posterior_bivariatezinba",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(marginal.mu), # double* marginal_mu
    as.double(marginal.size), # double* marginal_size
    as.double(marginal.beta), # double* marginal_beta
    as.double(size.coef), # double* size_coef
    as.double(mu.coef), # double* mu_coef
    as.double(beta0), # double* beta0
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N)) # double* proba
  # POST will hold the results

  posterior = matrix(z[[12]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[10]], ncol=N), proba=z[[11]])
  attr(posterior, "params") = params
  return(posterior)
}


# test function for the zinba copula
dzinba.copula <- function(x, y, fit, log=FALSE) {
  O = rbind(x, y)
  cormat = cov2cor(fit$sigma)
  z = .C("R_dzinbacopula",
    as.double(O),
    as.integer(length(x)),
    as.double(fit$marginal.x["size"]),
    as.double(fit$marginal.x["mu"]),
    as.double(fit$marginal.x["beta"]),
    as.double(fit$marginal.y["size"]),
    as.double(fit$marginal.y["mu"]),
    as.double(fit$marginal.y["beta"]),
    as.double(sqrt(fit$sigma[1,1])),
    as.double(sqrt(fit$sigma[2,2])),
    as.double(cormat[1,2]),
    0, # as.double(fit$mu[1]), # 0 is the old default
    0, # as.double(fit$mu[2]), # 0 is the old default
    double(length=length(x)))
  return(z[[14]])
}

dmvzinba.copula <- function(x, fit, log=FALSE) {
  O = t(x)
  z = .C("R_dmvzinbacopula",
    as.double(O),
    as.integer(ncol(O)),
    as.integer(nrow(O)),
    as.double(fit$marginals[,"size"]),
    as.double(fit$marginals[,"mu"]),
    as.double(fit$marginals[,"beta"]),
    as.double(fit$sigma),
    double(length=ncol(O)))
  return(z[[8]])
}

dmvn <- function(x, mu, sigma, log=FALSE) {
  O = t(x)
  z = .C("R_dmvn",
    as.double(O),
    as.integer(ncol(O)),
    as.integer(nrow(O)),
    as.double(mu),
    as.double(sigma),
    double(length=ncol(O)))
  return(z[[6]])
}


zinbacopula.hmm.posterior <- function(x, fits, max.it=500, eps=1e-4) {

  # number of states
  N = length(fits)
  size.x = sapply(fits, function(x) x$marginal.x["size"])
  mu.x = sapply(fits, function(x) x$marginal.x["mu"])
  beta.x = sapply(fits, function(x) x$marginal.x["beta"])
  size.y = sapply(fits, function(x) x$marginal.y["size"])
  mu.y = sapply(fits, function(x) x$marginal.y["mu"])
  beta.y = sapply(fits, function(x) x$marginal.y["beta"])
  sigma.x = sapply(fits, function(x) sqrt(x$sigma[1,1]))
  sigma.y = sapply(fits, function(x) sqrt(x$sigma[2,2]))
  rho = sapply(fits, function(x) cov2cor(x$sigma)[1,2])
  
  
  # length of observation sequence
  T = ncol(x) 
 
  #  
  z = .C("R_hmm_posterior_zinbacopula",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(size.x),
    as.double(mu.x),
    as.double(beta.x),
    as.double(size.y),
    as.double(mu.y),
    as.double(beta.y),
    as.double(sigma.x),
    as.double(sigma.y),
    as.double(rho),
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N)) # double* proba
  # POST will hold the results

  posterior = matrix(z[[15]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[16]], ncol=N), proba=z[[17]])
  attr(posterior, "params") = params
  return(posterior)
}

mvzinbacopula.hmm.posterior <- function(x, fits, max.it=500, eps=1e-4) {

  # number of states
  N = length(fits)
  size = sapply(fits, function(x) x$marginals[,"size"])
  mu = sapply(fits, function(x) x$marginals[,"mu"])
  beta = sapply(fits, function(x) x$marginal[,"beta"])
  sigma = sapply(fits, function(x) as.numeric(x$sigma))
  
  
  # length of observation sequence
  T = ncol(x) 
 
  #  
  z = .C("R_hmm_posterior_mvzinbacopula",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.integer(nrow(x)),
    as.double(size),
    as.double(mu),
    as.double(beta),
    as.double(sigma),
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N)) # double* proba
  # POST will hold the results

  posterior = matrix(z[[11]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[12]], ncol=N), proba=z[[13]])
  attr(posterior, "params") = params
  return(posterior)
}


hmm.baumwelch <- function(x, meanlog, sdlog, max.it=500, eps=1e-4, rho=0) {

  # number of states
  N = length(meanlog)
  stopifnot(length(sdlog) == N)
  # length of observation sequence
  T = length(x)
  
  # we also allow bivariate input if x is a 2-row matrix
  bivariate = is.matrix(x) && nrow(x) == 2
  if (bivariate) {
    # make sure that meanlog and sdlog are also 2-row matrices
    stopifnot(is.matrix(meanlog), nrow(meanlog) == 2,
              is.matrix(sdlog), nrow(sdlog) == 2)
    N = ncol(meanlog)
    T = ncol(x)
    if (length(rho) == 1) {
      rho = rep(rho, N)
    }
    stopifnot(length(rho) == N)
  }
  
  # void R_hmm_posterior(double* O, int* T, int* N, double* mu, double* sigma, int* iterationMAX, double* eps, double* post)
  z = .C("R_hmm_baumwelch",
    as.double(x + 1), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(meanlog), # double* mu
    as.double(sdlog), # double* sigma
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    as.integer(bivariate), # int* bivariate
    double(length=N*N), # double* A
    double(length=N), # double* proba
    as.double(rho)) # double* rho
  # POST will hold the results

  posterior = matrix(z[[8]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[10]], ncol=N), proba=z[[11]], meanlog=z[[4]], sdlog=z[[5]])
  attr(posterior, "params") = params
  return(posterior)
}


hmm.baumwelch.negbinom <- function(x, mu, size, max.it=500, eps=1e-4) {

  # number of states
  N = length(mu)
  stopifnot(length(size) == N)
  # length of observation sequence
  T = length(x)
  
  # void R_hmm_baumwelch_negbinom(double* O, int* T, int* N, double* mu, double* size, int* iterationMAX, double* eps, double* post, double* A, double* proba)

  z = .C("R_hmm_baumwelch_negbinom",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.double(mu), # double* mu
    as.double(size), # double* sigma
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N)) # double* proba
  # POST will hold the results

  posterior = matrix(z[[8]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[9]], ncol=N), proba=z[[10]], mu=z[[4]], size=z[[5]])
  attr(posterior, "params") = params
  return(posterior)
}


callRegions <- function(avg.posterior, cutoff, posterior.col="avg.posterior", expr.col="avg.expression") {
  require(GenomicRanges)

  called.regions = avg.posterior[which(avg.posterior[,posterior.col] > cutoff),]
  called.granges = GRanges(seqnames=called.regions$chrom, ranges=IRanges(start=called.regions$start, end=called.regions$end))
  elementMetadata(called.granges) = DataFrame(called.regions[,c(posterior.col, expr.col)])
  colnames(elementMetadata(called.granges)) = c(posterior.col, expr.col)
  
  # merge all neighboring regions
  reduced.granges = reduce(called.granges)
  # map bins to merged regions
  # called2reduced = as.matrix(findOverlaps(called.granges, reduced.granges))
  called2reduced = findOverlaps(called.granges, reduced.granges)
  called2reduced = cbind(queryHits(called2reduced), subjectHits(called2reduced))
  
  # compute global avg posterior
  reduced.posterior = tapply(elementMetadata(called.granges)[,posterior.col],
    called2reduced[,2], mean, na.rm=T)
  # compute global avg expression
  if (!is.null(expr.col)) {
    reduced.expression = tapply(elementMetadata(called.granges)[,expr.col],
      called2reduced[,2], mean, na.rm=T)
  } else {
    reduced.expression = rep(NA, length(reduced.granges))
  }
  elementMetadata(reduced.granges) = DataFrame(avg.posterior=reduced.posterior,
                   avg.expression=reduced.expression)
  return(reduced.granges)
}
