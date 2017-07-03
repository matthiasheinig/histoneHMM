zinba.mcmc <- function(x, N=10000) {
   z = .C("R_zinba_mcmc",
     as.double(x), # double* O,
     as.integer(length(x)), # int* N,
     as.double(0), # double* mu,
     as.double(0), # double* size,
     as.double(0), # double* beta,
     as.integer(N)) # int* iterations)
   return(c(mu=z[[3]], size=z[[4]], beta=z[[5]]))
}

logsum <- function(x) {
  m = max(x)
  return(log(sum(exp(x - m))) + m)
}

dzinba <- function(x, size, mu, beta, log=FALSE) {
  # beta is the inflation parameter
  if (beta < 0 || beta > 1) {
    return(rep(NA, length(x)))
  }
  alpha = 1 - beta
  is.zero = (x == 0)
  regular.density = dnbinom(x, size=size, mu=mu, log=log)
  if (log == FALSE) {
    density = alpha * regular.density
    density[is.zero] = density[is.zero] + beta
  } else {
    density = log(alpha) + regular.density
    density[is.zero] = sapply(density[is.zero], function(x) logsum(x, log(beta)))
  }
  return(density)
}

.pzinba <- function(x, size, mu, beta, lower.tail=TRUE) {
  # for values < 0 it returns 0
  cdf = stepfun(0:max(c(x, 1)), c(0, cumsum(dzinba(0:max(c(x, 1)), size, mu, beta))))
  p = cdf(x)
  if (!lower.tail) {
    p = 1 - p
  }
  return(p)
}

# there is also an analytical form of the cdf
pzinba <- function(x, size, mu, beta, lower.tail=TRUE) {
  # change the parameterization
  p = size / (size + mu)
  ## this is from mathworld.wolfram.com
  ## pbeta is actually the regularized incomplete beta function (see ?pbeta)
  cdf = pbeta(p, size, x + 1) 
  ## now we need to add the zero inflation part
  cdf = beta + (1 - beta) * cdf
  return(cdf)
}


rzinba <- function(n, size, mu, beta) {
  x = rnbinom(n, size=size, mu=mu)
  i = sample(n, round(beta * n))
  x[i] = 0
  return(x)
}

qzinba <- function(p, size, mu, beta, lower.tail=TRUE) {
  x = 0:5000
  cdf = cumsum(dzinba(x[-length(x)], size, mu, beta))
  inv.cdf = stepfun(cdf, x)
  q = inv.cdf(p)
  if (!lower.tail) {
    q = 1 - q
  }
  return(q)
}

fitzinba <- function(x, weight=1, start=NULL, fix.beta=FALSE) {
  # numerical fitting, we allow a weight to enable a 'soft' EM
  if (length(weight) == 1) {
    weight = rep(weight, length(x))
  }
  n = sum(weight)
 
  if (is.null(start)) {
    if (!fix.beta) {
      beta = sum(weight[x == 0], na.rm=T) / n + 0.000001 # make sure it is not 0
    } else {
      beta = 0
    }
    # what follows is modified from fitdistr (MASS)
    m = x %*% weight / n
    v = sum(weight * (x - m)^2) / n
    size = ifelse(v > m, m^2/(v - m), 100)
    # we just add the beta
    start = list(size=size, mu=m, beta=beta)
    cat("initialize:")
    print(start)
    # objective function is the -loglikelihood (optim minimizes)
    # actually this is the expected loglikelihood (using the weights)
  } else {
    start = as.list(start)
  }
  if (fix.beta) {
    beta = start[["beta"]]
    start = start[-match("beta", names(start))]
    obj = function(par) -sum(weight * log(dzinba(x, par[["size"]], par[["mu"]], beta))) / n
  } else {
    obj = function(par) -sum(weight * log(dzinba(x, par[["size"]], par[["mu"]], par[["beta"]]))) / n
  }
  fit = optim(start, obj)
  if (fit$convergence != 0) {
    warning("fitzinba did not converge")
  }
  return(fit$par)
}

.fitzinba <- function(x, weight=1, start=NULL, fix.beta=FALSE) {
  if (length(weight) == 1) {
    weight = rep(weight, length(x))
  }
  if (!is.null(start)) {
    mu = start["mu"]
    size = start["size"]
    beta = start["beta"]
  } else {
    mu = 0
    size = 0
    beta = 0
  }
  if (fix.beta) {
    ## actually fit a negative binomial instead
    z = .C("R_fitnegbinom", as.double(x), as.double(weight), as.integer(length(x)), as.double(mu), as.double(size))
    fit = c(mu=z[[4]], size=z[[5]], beta=0)
    return(fit)
  } else {
    # call to the c function void R_fitzinba(double* x, double* w, int* n, double* mu, double* size, double* beta)
    z = .C("R_fitzinba", as.double(x), as.double(weight), as.integer(length(x)), as.double(mu), as.double(size), as.double(beta))
    fit = c(mu=z[[4]], size=z[[5]], beta=z[[6]])
    return(fit)
  }
}


if (F) {
  # test this
  library(histoneHMM)
  mu = 5
  size = 10
  beta = 0.1
  estimates = t(sapply(1:100, function(i) {
    x = rzinba(100, mu=mu, size=size, beta=beta)
    return(.fitzinba(x))
  }))

  # also test the em
  x = c(rzinba(100, mu=mu, size=size, beta=beta), rzinba(100, mu=mu+10, size=size, beta=beta))


  x = rnbinom(1000, mu=5, size=10)
  .fitzinba(x, fix.beta=T)

  partial <- function(y, size) {
    N = length(y)
    p = N * size / (N * size + sum(y))
    return(N * log(p) - N * digamma(size) + sum(digamma(size + y)))
  }
  partial2 <- function(y, w, size) {
    sum.w = sum(w)
    sum.wy = sum(y * w)
    p = size * sum.w / (sum.wy + size * sum.w)
    l = sum(w * (digamma(y + size) - digamma(size) + log(p)))
    return(l)
  }

}


zinba.qqplot <- function(x, size, mu, beta) {
  # no idea how to compute the quantiles more efficiently..
  qfun = stepfun(cumsum(dzinba(0:max(x), size, mu, beta)), 0:(max(x)+1))
  sx = sort(x)
  sy = qfun(ecdf(sx)(sx))
  plot(sx, sy, xlab="observed", ylab="theoretical")
  abline(a=0, b=1)
}



logdensityplot <- function(x, size, mu, beta) {
  sx = sort(x)
  plot(x=sx, y=1-ecdf(x)(sx), log="xy", type="l", xlab="x", ylab="1 - CDF(x)")
  lines(x=0:max(x), 1 - pzinba(0:max(x), size, mu, beta), col="red")
  legend("topright", lty=1, col=c("black", "red"), legend=c("data", "fitted"))
}

# since we have a discrete distribution we need a modified version of the
# probability integral transformation (suggested also in)
# The Journal of Risk Model Validation (3–10) Volume 3/Number 2, Summer 2009
# A note on the Berkowitz test with discrete distributions
# Alfred Hamerle and Kilian Plank

zinba.pit <- function(x, size, mu, beta, lower.tail=TRUE) {
  # first just the regular cdf
  cdf.val = cumsum(dzinba(0:max(c(x, 1)), size, mu, beta))
  cdf = stepfun(1:max(c(x, 1)), cdf.val)
  p = cdf(x)
  if (F) {
  # now we need to distribute the tied observations uniformly in the intervals
  # defined by cdf.val
  intervals = c(0, cdf.val, 1)
  # sometimes it happens that the cdf value for the max is numerically 1
  # this causes trouble with the findInterval function, so we remove dups
  intervals = intervals[!duplicated(intervals)]
  i = findInterval(p, intervals[-1])
  p = runif(length(p), min=intervals[i], max=intervals[i+1])
  }
  
  if (!lower.tail) {
    p = 1 - p
  }
  return(p)
}


fit.zinba.copula <- function(x, y, weight=1, marginal.x=NULL, marginal.y=NULL) {
  if (is.null(marginal.x)) {
    marginal.x = fitzinba(x, weight)
  }
  if (is.null(marginal.y)) {
    marginal.y = fitzinba(y, weight)
  }
  
  if (length(weight) == 1) {
    weight = rep(weight, length(x))
  }

  # tranform to p-values [0, 1]
  p.x = zinba.pit(x, marginal.x["size"], marginal.x["mu"], marginal.x["beta"])
  p.y = zinba.pit(y, marginal.y["size"], marginal.y["mu"], marginal.y["beta"])

  # transform to z-values
  z.x = qnorm(p.x)
  z.y = qnorm(p.y)

  # get the covariance matrix (we only use the non-zero observations)
  use = !(x == 0 | y == 0) & !(is.na(z.x) | is.na(z.y)) & (z.x < Inf & z.y < Inf)
  z.mat = cbind(z.x[use], z.y[use])
  mu = colMeans(z.mat)
  z.mat = (z.mat - rep(mu, each=nrow(z.mat))) * sqrt(weight[use])
  n = sum(weight[use])
  sigma = 1 / n * t(z.mat) %*% z.mat

  return(list(marginal.x=marginal.x, marginal.y=marginal.y, sigma=sigma, mu=mu))
}


fit.zinba.copula2 <- function(x, weight=1, marginals=NULL) {
  if (is.null(marginals)) {
    marginals = t(apply(x, 2, fitzinba, weight))
  }
  
  if (length(weight) == 1) {
    weight = rep(weight, nrow(x))
  }
  n = sum(weight)

  # tranform to p-values [0, 1]
  p = sapply(1:ncol(x), function(i) zinba.pit(x[,i], marginals[i,"size"], marginals[i,"mu"], marginals[i,"beta"]))
 
  # transform to z-values
  z.mat = apply(p, 2, qnorm)

  # get the covariance matrix
  mu = colMeans(z.mat)
  z.mat = (z.mat - rep(mu, each=nrow(z.mat))) * sqrt(weight)
  sigma = 1 / n * t(z.mat) %*% z.mat

  return(list(marginals=marginals, sigma=sigma, mu=mu))
}



# this version computes the density from the CDF:
# P(X=x, Y=y) = P(X<x, Y<y) - P(X<x-1, Y<y) - P(X<x, Y<y-1) + P(X<x-1, Y<y-1)
.dzinba.copula <- function(x, y, fit, log=FALSE, algorithm=NULL) {
  require(mvtnorm)  

  # get rid of duplicated pairs to avoid expensive mv normal CDF calls
  # pair = paste(x, y)
  # o = order(pair)
  # pair = pair[o]
  # uniq = unique(pair)
  # p2u = match(pair, uniq)
  # uniq = unique(p2u)
  # 
  # x = x[o][uniq]
  # y = y[o][uniq]
  
  # tranform to p-values
  upper.p.x = pzinba(x, fit$marginal.x["size"], fit$marginal.x["mu"], fit$marginal.x["beta"])
  upper.p.y = pzinba(y, fit$marginal.y["size"], fit$marginal.y["mu"], fit$marginal.y["beta"])
 
  # transform to z-values
  upper.z.x = qnorm(upper.p.x)
  upper.z.y = qnorm(upper.p.y)

  # tranform to p-values
  lower.p.x = pzinba(x - 1, fit$marginal.x["size"], fit$marginal.x["mu"], fit$marginal.x["beta"])
  lower.p.y = pzinba(y - 1, fit$marginal.y["size"], fit$marginal.y["mu"], fit$marginal.y["beta"])
 
  # transform to z-values
  lower.z.x = qnorm(lower.p.x)
  lower.z.y = qnorm(lower.p.y)
 
  # compute the CDF in intervals (x-1, y-1) to (x, y)

  # this is a monte carlo based integration method
  if (is.null(algorithm)) {
    algorithm = GenzBretz()
  }
  d = apply(cbind(lower.z.x, lower.z.y, upper.z.x, upper.z.y), 1, function(x) pmvnorm(lower=x[1:2], upper=x[3:4], mean=rep(0, 2), sigma=fit$sigma, algorithm=algorithm))
  
  if (log == TRUE) {
    d = log(d)
  }

  # get back to the original order
  # dd = d[p2u]
  # dd[o] = dd
  dd = d
  
  return(dd)
}




pzinba.copula <- function(x, y, fit) {
  require(mvtnorm)
  
  # tranform to p-values
  p.x = pzinba(x, fit$marginal.x["size"], fit$marginal.x["mu"], fit$marginal.x["beta"])
  p.y = pzinba(y, fit$marginal.y["size"], fit$marginal.y["mu"], fit$marginal.y["beta"])
 
  # transform to z-values
  z.x = qnorm(p.x)
  z.y = qnorm(p.y)
 
  # compute the CDF
  p = pmvnorm(lower=rep(-Inf, 2), upper=cbind(z.x, z.y), mean=rep(0, 2), sigma=fit$sigma)
  return(p)
}


rzinba.copula <- function(n, fit) {
  #TODO: maybe do lib = require(mvtnorm); if(lib == FALSE){install.package("mvtnorm"); require(mvtnorm);}
  require(mvtnorm)

  # generate observations in the gaussian space
  z = rmvnorm(n, mean=fit$mu, sigma=fit$sigma)

  # transfor to uniform variates [0, 1]
  p = apply(z, 2, pnorm)
  
  # transform back to the count space using the marginals
  if (!is.null(fit$marginal.x)) {
    x = qzinba(p[,1], fit$marginal.x["size"], fit$marginal.x["mu"], fit$marginal.x["beta"])
    y = qzinba(p[,2], fit$marginal.y["size"], fit$marginal.y["mu"], fit$marginal.y["beta"])
    return(cbind(x, y))
  } else {
    return(sapply(1:ncol(z), function(i) qzinba(p[,i], fit$marginals[i,"size"], fit$marginals[i,"mu"], fit$marginals[i,"beta"])))
  }
}


# little helper to save and load zinbacopula fits
fit.to.col <- function(fit) {
  res = c(fit$marginal.x, fit$marginal.y, fit$sigma[upper.tri(fit$sigma, diag=T)], fit$mu)
  names(res) = c(paste("marginal.x", names(fit$marginal.x), sep="."),
         paste("marginal.y", names(fit$marginal.x), sep="."),
         paste("sigma", c("11", "12", "22"), sep="."),
         paste("mu", 1:2, sep="."))
  return(res)
}

col.to.fit <- function(col) {
  sigma = matrix(nrow=2, ncol=2)
  sigma[upper.tri(sigma, diag=T)] = col[grep("sigma", names(col))]
  sigma[2,1] = sigma[1,2]
  res = list(marginal.x=col[grep("marginal.x", names(col))],
    marginal.y=col[grep("marginal.y", names(col))],
    sigma=sigma, mu=col[grep("^mu", names(col))])
  names(res$marginal.x) = gsub("marginal.x.", "", names(res$marginal.x))
  names(res$marginal.y) = gsub("marginal.y.", "", names(res$marginal.y))
  return(res)
}

