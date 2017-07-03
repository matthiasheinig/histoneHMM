simulate.lognormal <- function(T, A, proba, meanlog, sdlog, rho=0) {
  p = runif(T)
  proba = cumsum(proba)
  A = t(apply(A, 1, cumsum))
  # initial state
  states = numeric(length=T)
  states[1] = findInterval(p[1], proba, rightmost.closed=T) + 1
  # all other states
  sapply(2:T, function(i) {
    states[i] <<- findInterval(p[i], A[states[i - 1],], rightmost.closed=T) + 1
  })
  # observations from the states
  if (is.matrix(meanlog) && nrow(meanlog) == 2 && is.matrix(sdlog) && nrow(sdlog) == 2) {
    require(tmvtnorm) # truncated multivariate normal (>0!!)
    # in the bivariate case we draw from a multivariate normal and then exp it
    O = matrix(ncol=2, nrow=T)
    for (state in 1:nrow(A)) {
      idx = states == state
      covariance = rho[state] * prod(sdlog[,state])
      sigma = matrix(covariance, nrow=2, ncol=2)
      diag(sigma) = sdlog[,state]^2
      mvnorm = rtmvnorm(sum(idx), mean=meanlog[,state], sigma=sigma, lower=rep(0, 2))
      O[idx,] = exp(mvnorm)
    }
  } else {
    O = rlnorm(T, meanlog=meanlog[states], sdlog=sdlog[states])
  }
  return(cbind(states, O))
}


simulate.zinba.copula <- function(T, A, proba, fits) {
  p = runif(T)
  proba = cumsum(proba)
  A = t(apply(A, 1, cumsum))
  # find out how many dimensions we have
  pp = 2
  if (!is.null(fits[[1]]$marginals)) {
    pp = nrow(fits[[1]]$marginals)
  }
  # initial state
  states = numeric(length=T)
  states[1] = findInterval(p[1], proba, rightmost.closed=T) + 1
  # all other states
  sapply(2:T, function(i) {
    states[i] <<- findInterval(p[i], A[states[i - 1],], rightmost.closed=T) + 1
  })
  # observations from the states
  O = matrix(ncol=pp, nrow=T)
  for (state in 1:nrow(A)) {
    idx = states == state
    obs = rzinba.copula(sum(idx), fits[[state]])
    O[idx,] = obs
  }
  return(cbind(states, O))
}


## this should be used and developped further
simulate.HMM <- function (T, A, proba, models) {
    p = runif(T)
    proba = cumsum(proba)
    A = t(apply(A, 1, cumsum))
    states = numeric(length = T)
    states[1] = findInterval(p[1], proba, rightmost.closed = T) + 1
    sapply(2:T, function(i) {
        states[i] <<- findInterval(p[i], A[states[i - 1], ], 
            rightmost.closed = T) + 1
    })

    O = matrix(ncol = 1, nrow = T)
    for (state in 1:nrow(A)) {
        idx = states == state
        O[idx, ] = drawFrom(models[[state]], sum(idx))
    }
    return(cbind(states, O))
}
