# for a comparison between using expression data and just using the signal
# we need an EM algorithm to estimate a log normal mixture
library(methods)

# helper for more numerical precision when summing up log likelihoods
logsum <- function(x) {
  m = max(x)
  return(log(sum(exp(x - m))) + m)
}

# distributions are defined in a object oriented way
# they have a logLikelihood and a fit function

# dummy generic logLikelihood function (returns -Inf)
logLikelihood <- function(model, x) return (-Inf)
setGeneric("logLikelihood")

# dummy generic fit function (just returns the model)
fit <- function(model, x, logExpectation) {
  return(model)
}
setGeneric("fit")

drawFrom <- function(model, N) {
    print("drawFrom dummy")
    return(rep(NA, N))
}
setGeneric("drawFrom")

CDF <- function(model, x) {
  print("CDF dummy")
  return(rep(NA, length(x)))
}
setGeneric("CDF")

## defualt is just to return the total length of all slots of the object
nparam <- function(model) {
  slots <- slotNames(model)
  return(sum(sapply(slots, function(sname) length(slot(model, sname)))))
}
setGeneric("nparam")

# the superclass is Distribution
Distribution <- setClass("Distribution")

# allow plotting of the densities and mixture
plotDensity <- function(model, x=NA, xlim=NULL, add=FALSE, col="black", alpha=1) {
  print("plotDensity (generic)")
}
setGeneric("plotDensity")

plotDensity.Distribution <- function(model, x=NULL, xlim=NULL, add=FALSE, col="black", alpha=1) {
  # print("plotDensity.Distribution")
  if (is.null(x)) {
    if (is.null(xlim))
      stop("please specify either data or xlim")
  } else {
    if (is.null(xlim)) {
      xlim = range(x)
    } else {
      x = x[xlim[1] <= x & x <= xlim[2]]
    }
  }
  xx = seq(xlim[1], xlim[2])
  d = exp(logLikelihood(model, xx)) * alpha
  if (!is.na(x)[1] && !add) {
    hist(x, freq=F, add=add, breaks=xlim[2], main="", xlab="", ylab="")
  }
  if (!add || !is.na(x)) {
    lines(xx, d, col=col, lwd=2)
  } else {
    plot(xx, d, type="l", col=col, lwd=2, main="", xlab="", ylab="")
  }
}
setMethod("plotDensity", signature=c(model="Distribution", x="numeric", xlim="numeric", add="logical", col="character", alpha="numeric"), plotDensity.Distribution)





# for now we just put the lognormal distribution
LogNormal <- setClass("LogNormal", representation=list(meanlog="numeric", sdlog="numeric"))

logLikelihood.LogNormal <- function(model, x) {
  return(dlnorm(x, meanlog=model@meanlog, sdlog=model@sdlog, log=T))
}
setMethod("logLikelihood", signature=c(model="LogNormal", x="numeric"), logLikelihood.LogNormal)

fit.LogNormal <- function(model, x, logExpectation) {
  # adapted code from MASS::fitdistr
  # we need to add a weighting by the Expectation
  expect = exp(logExpectation)
  lx = log(x)
  n = sum(expect)
  meanlog = sum(lx * expect) / n
  sdlog = sqrt(sum(expect * (lx - meanlog)^2) / n)
  newmodel = new("LogNormal")
  newmodel@meanlog = meanlog
  newmodel@sdlog = sdlog
  return(newmodel)
}
setMethod("fit", signature=c(model="LogNormal", x="numeric", logExpectation="numeric"), fit.LogNormal)

# we also define this for a zero inflated negative binomial
Zinba <- setClass("Zinba", representation=list(size="numeric", mu="numeric", beta="numeric"), contains="Distribution")

logLikelihood.Zinba <- function(model, x) {
  return(log(dzinba(x, size=model@size, mu=model@mu, beta=model@beta)))
}
setMethod("logLikelihood", signature=c(model="Zinba", x="numeric"), logLikelihood.Zinba)

fit.Zinba <- function(model, x, logExpectation) {
  if (length(model@size) == 0) {
    start = NULL
  } else {
    start = c(model@size, model@mu, model@beta)
    names(start) = c("size", "mu", "beta")
  }
  par = fitzinba(x, weight=exp(logExpectation), start=start)
  # newmodel = Zinba()
  newmodel = new("Zinba")
  newmodel@size = par["size"]
  newmodel@mu = par["mu"]
  newmodel@beta = par["beta"]
  return(newmodel)
}
setMethod("fit", signature=c(model="Zinba", x="numeric", logExpectation="numeric"), fit.Zinba)

drawFrom.Zinba <- function(model, N) {
    return(rzinba(N, model@size, model@mu, model@beta))
}
setMethod("drawFrom", signature=c(model="Zinba", N="numeric"), drawFrom.Zinba)

Zero <- setClass("Zero", representation=list(p="numeric"), contains="Distribution")

fit.Zero <- function(model, x, logExpectation) {
    return(model)
}
setMethod("fit", signature=c(model="Zero", x="numeric", logExpectation="numeric"), fit.Zero)

logLikelihood.Zero <- function(model, x) {
  ll = rep(-Inf, length(x))
  ll[x == 0] = 0
  return(ll)
}
setMethod("logLikelihood", signature=c(model="Zero", x="numeric"), logLikelihood.Zero)

drawFrom.Zero <- function(model, N) {
    return(rep(0, N))
}
setMethod("drawFrom", signature=c(model="Zero", N="numeric"), drawFrom.Zero)

nparam.Zero <- function(model) return(0)
setMethod("nparam", signature=c(model="Zero"), nparam.Zero)


# we also define this for a zero inflated negative binomial using a C implementation
CZinba <- setClass("CZinba", representation=list(size="numeric", mu="numeric", beta="numeric"), contains="Distribution")

logLikelihood.CZinba <- function(model, x) {
  return(log(dzinba(x, size=model@size, mu=model@mu, beta=model@beta)))
}
setMethod("logLikelihood", signature=c(model="CZinba", x="numeric"), logLikelihood.CZinba)

fit.CZinba <- function(model, x, logExpectation) {
  if (length(model@size) == 0) {
    start = NULL
  } else {
    start = c(model@size, model@mu, model@beta)
    names(start) = c("size", "mu", "beta")
  }
  par = .fitzinba(x, weight=exp(logExpectation), start=start)
  # newmodel = CZinba()
  newmodel = new("CZinba")
  newmodel@size = par["size"]
  newmodel@mu = par["mu"]
  newmodel@beta = par["beta"]
  return(newmodel)
}
setMethod("fit", signature=c(model="CZinba", x="numeric", logExpectation="numeric"), fit.CZinba)




# just the regular negative binomial (special case of the zero inflated)
Negbinom <- setClass("Negbinom", representation=list(size="numeric", mu="numeric"), contains="Distribution")

logLikelihood.Negbinom <- function(model, x) {
  return(dnbinom(x, size=model@size, mu=model@mu, log=T))
}
setMethod("logLikelihood", signature=c(model="Negbinom", x="numeric"), logLikelihood.Negbinom)

## fit.Negbinom <- function(model, x, logExpectation) {
##   if (length(model@size) == 0) {
##     start = NULL
##   } else {
##     start = c(model@size, model@mu, 0)
##     names(start) = c("size", "mu", "beta")
##   }
##   par = fitzinba(x, weight=exp(logExpectation), start, fix.beta=TRUE)
##   # newmodel = Zinba()
##   newmodel = new("Negbinom")
##   newmodel@size = par["size"]
##   newmodel@mu = par["mu"]
##   return(newmodel)
## }

fit.Negbinom <- function(model, x, logExpectation) {
  par = .fitzinba(x, weight=exp(logExpectation), fix.beta=T)
  # newmodel = Zinba()
  newmodel = new("Negbinom")
  newmodel@size = par["size"]
  newmodel@mu = par["mu"]
  return(newmodel)
}
setMethod("fit", signature=c(model="Negbinom", x="numeric", logExpectation="numeric"), fit.Negbinom)

drawFrom.Negbinom <- function(model, N) {
    return(rnbinom(N, model@size, model@mu))
}
setMethod("drawFrom", signature=c(model="Negbinom", N="numeric"), drawFrom.Negbinom)

CDF.Negbinom <- function(model, x) {
  return(pnbinom(x, size=model@size, mu=model@mu))
}
setMethod("CDF", signature=c(model="Negbinom", x="numeric"), CDF.Negbinom)


## Gamma distribution
Gamma <- setClass("Gamma", representation=list(shape="numeric", rate="numeric"), contains="Distribution")

logLikelihood.Gamma <- function(model, x) {
  return(dgamma(x, shape=model@shape, rate=model@rate, log=T))
}
setMethod("logLikelihood", signature=c(model="Gamma", x="numeric"), logLikelihood.Gamma)

fit.Gamma <- function(model, x, logExpectation) {
  ## this is modified from MASS::fidistr
  if (any(x <= 0)) 
    warning("gamma values must be > 0")
  
  weight <- exp(logExpectation)
  m <- x %*% weight / n
  v <- sum(weight * (x - m)^2) / n

  start <- list(shape = m^2/v, rate = m/v)
  control <- list(parscale = c(1, start$rate))

  ## define the objective function
  obj <- function(par) -sum(weight * dgamma(x[x > 0], shape=par[["shape"]], rate=par[["rate"]], log=TRUE)) / n

  fit <- optim(start, obj)
  if (fit$convergence != 0) {
    warning("fitzinba did not converge")
  }
  par <- fit$par

  newmodel = new("Gamma")
  newmodel@shape = par["shape"]
  newmodel@rate = par["rate"]
  return(newmodel)
}
setMethod("fit", signature=c(model="Gamma", x="numeric", logExpectation="numeric"), fit.Gamma)

drawFrom.Gamma <- function(model, N) {
  return(rgamma(N, shape=model@shape, rate=model@rate))
}
setMethod("drawFrom", signature=c(model="Gamma", N="numeric"), drawFrom.Gamma)

CDF.Gamma <- function(model, x) {
  return(pnbinom(x, shape=model@shape, rate=model@rate))
}
setMethod("CDF", signature=c(model="Gamma", x="numeric"), CDF.Gamma)


## the EM behaves strangely somtimes, so we try truncated distributions
## beta is the normalization factor that we need to get a real density
## it is 1 / cdf(x) for left=F and 1 / (1-cdf(x)) for left=T
setClass("TruncatedDistribution", representation=list(distribution="Distribution", cutpoint="numeric", left="logical", beta="numeric"), contains="Distribution")
TruncatedDistribution <- function(distribution, cutpoint, left) {
  beta = CDF(distribution, cutpoint)
  if (left) {
    beta = 1 - beta
  } 
  beta = 1 / beta
  ## some problems with this, so we just use 1 for the moment
  beta = 1
  td = new("TruncatedDistribution", distribution=distribution, cutpoint=cutpoint, left=left, beta=beta)
  return(td)
}


truncate <- function(model, x) {
  if (model@left) {
    truncate = x <= model@cutpoint
  } else {
    truncate = x > model@cutpoint
  }
  return(truncate)
}
# setMethod("truncate", signature=c(model="TruncatedDistribution", x="numeric"))

logLikelihood.TruncatedDistribution <- function(model, x) {
  ll = logLikelihood(model@distribution, x)
  ll = ll + log(model@beta)
  truncated = truncate(model, x)
  ll[truncated] = -Inf
  return(ll)
}
setMethod("logLikelihood", signature=c(model="TruncatedDistribution", x="numeric"), logLikelihood.TruncatedDistribution)

fit.TruncatedDistribution <- function(model, x, logExpectation) {
  truncated = truncate(model, x)
  lexp = logExpectation
  lexp[truncated] = -Inf
  
  newdist = fit(model@distribution, x, lexp)
  newmodel = TruncatedDistribution(newdist, model@cutpoint, model@left)

  return(newmodel)
}
setMethod("fit", signature=c(model="TruncatedDistribution", x="numeric", logExpectation="numeric"), fit.TruncatedDistribution)

drawFrom.TruncatedDistribution <- function(model, N) {
    return(rnbinom(N, model@size, model@mu))
}
setMethod("drawFrom", signature=c(model="TruncatedDistribution", N="numeric"), drawFrom.TruncatedDistribution)





# copula based bivariate negative binomial
ZinbaCopula <- setClass("ZinbaCopula", representation=list(marginal.x="numeric", marginal.y="numeric", sigma="matrix"), contains="Distribution")

logLikelihood.ZinbaCopula <- function(model, x) {
  require(histoneHMM)
  lobj = list(marginal.x=model@marginal.x, marginal.y=model@marginal.y, sigma=model@sigma)
  return(log(histoneHMM::dzinba.copula(x[,1], x[,2], lobj)))
}
setMethod("logLikelihood", signature=c(model="ZinbaCopula", x="matrix"), logLikelihood.ZinbaCopula)

fit.ZinbaCopula <- function(model, x, logExpectation) {
  par = fit.zinba.copula(x[,1], x[,2], weight=exp(logExpectation))
  newmodel = new(ZinbaCopula)
  newmodel@marginal.x = par$marginal.x
  newmodel@marginal.y = par$marginal.y
  newmodel@sigma = par$sigma
  return(newmodel)
}
setMethod("fit", signature=c(model="ZinbaCopula", x="matrix", logExpectation="numeric"), fit.ZinbaCopula)



EMResult <- setClass("EMResult", representation=list(models="list", logLikelihood="matrix", dataLogLikelihood="numeric", proportions="numeric"), contains="Distribution")

plotDensity.EMResult <- function(model, x=NULL, xlim=NULL, add=FALSE, col="black", alpha=1) {

  # plot the components and compute the mixture density
  xx = seq(xlim[1], xlim[2])
  d = rep(0, length(xx))
  for (i in 1:length(model@models)) {
    plotDensity(model@models[[i]], x, xlim, add=((i!=1) || (i==1 && add)), col=palette()[i+1], alpha=model@proportions[i])
    d = d + exp(logLikelihood(model@models[[i]], xx)) * model@proportions[i]
  }

  if (!add || !is.null(x)) {
    lines(xx, d, col=col, lwd=2)
  } else {
    plot(xx, d, type="l", col=col, lwd=2)
  }
  if (!is.null(names(model@models))) {
    legend("right", lty=1, lwd=2, col=c(col, palette()[1 + 1:length(model@models)]), legend=c("mixture", names(model@models)))
  }
}
setMethod("plotDensity", signature=c(model="EMResult", x="numeric", xlim="numeric", add="logical", col="character", alpha="numeric"), plotDensity.EMResult)

logLikelihood.EMResult <- function(model, x) {
  ll = sapply(model@models, logLikelihood, x)
  ll = ll + rep(log(model@proportions), each=length(x))
  ll = apply(ll, 1, logsum)
  return(ll)
}
setMethod("logLikelihood", signature=c(model="EMResult", x="numeric"), logLikelihood.EMResult)

nparam.EMResult <- function(model) {
  nparam.submodels <- sapply(model@models, nparam)
  ## proportions sum up to one so we actually have nr of components - 1 proportion parameters
  return(sum(nparam.submodels) + length(model@models) - 1)
}
setMethod("nparam", signature=c(model="EMResult"), nparam.EMResult)

EMPosterior <- function(emfit, x) {
  posterior = sapply(emfit@models, logLikelihood, x) + rep(log(emfit@proportions), each=length(x))
  posterior = posterior - apply(posterior, 1, logsum)
  posterior = exp(posterior)
  return(posterior)
}


## compute AIC and BIC for the mixture models
AIC <- function(emfit) {
  return(2 * nparam(emfit) - 2 * emfit@dataLogLikelihood)
}

BIC <- function(emfit) {
  return(log(nrow(emfit@logLikelihood)) * nparam(emfit) - 2 * emfit@dataLogLikelihood)
}


# iteratively estimate a mixture
# also accept preinitialized models
em <- function(x, ncomp, prop=NULL, maxit=100, eps=1e-4, model.constructor="Zinba", models=NULL) {
  n = length(x)
  if (is.matrix(x) || is.data.frame(x)) {
    n = nrow(x)
  }
  if (is.null(prop)) {
    prop = rep(1 / ncomp, ncomp)
  }
  logProp = log(prop)

  if (length(model.constructor) == 1) {
    model.constructor = rep(model.constructor, ncomp)
  } else if (length(model.constructor) != ncomp) {
    stop("model.constructor has to have length 1 or ncomp")
  }
  
  # initialize the models
  if (is.null(models)) {
    models = list()
    # we just use quantiles to cut the data
    if (is.matrix(x) || is.data.frame(x)) {
      # for simplicity just cut along the first axis
      cutted = findInterval(x[,1], quantile(x[,1], prob=(0:ncomp)/ncomp),  rightmost.closed=T)
    } else {
      cutted = findInterval(x, quantile(x, prob=(0:ncomp)/ncomp),  rightmost.closed=T)
    }

    # browser()
    cat("Initializing..\n")
    for (m in 1:ncomp) {
      model = new(model.constructor[m])
      # make the weights a mixture of the quantile based assignment but also
      # distribute a bit of weight on the other classes
      alpha = 0.9
      weight = alpha * as.numeric(cutted == m)
      weight[weight == 0] = 1 - alpha
      model = fit(model, x, log(weight))
      models[[m]] = model
    }
  } else {
    cat("Using preinitialized models\n")
  }
  

  oldLogLik = NA
  # iterate E and M step
  for (i in 1:maxit) {
    #
    # E step
    #
    logExpectation = sapply(models, logLikelihood, x) # p_i(x)
    logExpectation = logExpectation + rep(logProp, each=n) # a_i p_i(x)
    logLikByElement = apply(logExpectation, 1, logsum) # p(x)
    logExpectation = logExpectation - logLikByElement # a_i p_i(x) / p(x)
    logLik = sum(logLikByElement)

    # check convergence using the delta in log likelihood
    cat("E-step: logLik =", logLik, "delta =", logLik - oldLogLik, "\n")
    if (!is.na(oldLogLik) && logLik - oldLogLik < eps) {
      if (logLik - oldLogLik < 0) {
        warning("Problem in the EM algorithm: likelihood is decreasing!")
      }
      break
    }
    oldLogLik = logLik

    #
    # M step
    #

    # new estimate for the proportions
    logProb = apply(logExpectation, 2, logsum) - log(n)
    
    # new estimates for the models
    for (m in 1:ncomp) {
      models[[m]] = fit(models[[m]], x, logExpectation[,m])
    }

    cat("M-step: proportions =", exp(logProb), "\nmodels:\n")
    print(models)
    
    res = new("EMResult")
    res@models = models
    res@logLikelihood = logExpectation
    res@dataLogLikelihood = logLik
    res@proportions = exp(logProb)
  }
  return(res)
}



if (F) {

  # test this stuff
  x = rnbinom(1000, size=10, mu=5)
  # inflate the zeros
  x = c(x, rep(0, 100))

  y = rnbinom(1000, size=10, mu=10)
  y = c(y, rep(0, 20))
  
  fitzinba(x)
  fitzinba(y)

  # ok this works pretty nicely

  w1 = c(rep(0.9, length(x)), rep(0.1, length(y)))
  w2 = c(rep(0.1, length(x)), rep(0.9, length(y)))
  fitzinba(c(x, y), w1)
  fitzinba(c(x, y), w2)

  
  # now check the em
  emfit = em(c(x, y), ncomp=2, model.constructor="Zinba")

  # test with real data
  load("tmp/diff-data.RData")

  k27 = c(high1, low1)
  emfit = em(k27[k27 < 500], ncomp=2, model.constructor="Zinba")

  H3K27me3 = k27[k27 < 500]
  count = 0:max(H3K27me3)
  d1 = exp(logLikelihood(emfit@models[[1]], count))
  d2 = exp(logLikelihood(emfit@models[[2]], count))

  dmix = d1 * emfit@proportions[1] + d2 * emfit@proportions[2]
 
  hist(H3K27me3, breaks=50, freq=F)
  lines(count, d1 * emfit@proportions[1], col="green", lwd=2)
  lines(count, d2 * emfit@proportions[2], col="blue", lwd=2)
  lines(count, dmix, col="black", lwd=2)


  # test the copula model
  library(mvtnorm)
  sigma = matrix(c(1, 0.6,
                   0.6, 1), byrow=T, ncol=2)
  z = rmvnorm(1000, sigma=sigma)
  z[,1] = pnorm(z[,1])
  z[,2] = pnorm(z[,2])
  
  x = z
  x[,1] = qnbinom(x[,1], size=10, mu=10)
  x[,2] = qnbinom(x[,2], size=10, mu=10)

  y = z
  y[,1] = qnbinom(y[,1], size=30, mu=20)
  y[,2] = qnbinom(y[,2], size=30, mu=20)

  data = rbind(x, y)
  label = factor(rep(c("x", "y"), each=1000))

  emfit = em(data, ncomp=2, model.constructor=ZinbaCopula)

  # simulate data from two clearly seperable classes
  sigma1 = matrix(c(1, 0.6,
                    0.6, 1), byrow=T, ncol=2)

  angle = 1/8 * 2* pi

  rotation = matrix(c(cos(angle), -sin(angle),
                      sin(angle), cos(angle)), byrow=T, ncol=2)

  sigma2 = t(rotation) %*% sigma %*% rotation 

  
  x = rzinba.copula(1000, fit=list(mu=c(0, 0), sigma=sigma2, marginal.x=c(size=10, mu=10, beta=0.1), marginal.y=c(size=10, mu=10, beta=0.1)))
  
  y = rzinba.copula(1000, fit=list(mu=c(0, 0), sigma=sigma, marginal.x=c(size=20, mu=30, beta=0.001), marginal.y=c(size=20, mu=30, beta=0.001)))

  data = rbind(x, y)

  plot(data, col=label)

  emfit = em(data, ncomp=2, model.constructor=ZinbaCopula)

  
  ## test the gamma model
  n <- 500
  x <- c(rgamma(n, shape=3, rate=2), rgamma(n, shape=5, rate=1))
  
  hist(x, breaks=50)
  
  emfit <- em(x, ncomp=2, model.constructor=c("Gamma", "Gamma"))
  
  ## test gamma with zero inflation
  y <- c(rep(0, 50), x)
  emfit <- em(y, ncomp=3, model.constructor=c("Zero", "Gamma", "Gamma"))
  
  BIC(emfit)
  AIC(emfit)
}
