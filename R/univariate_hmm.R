univariate.hmm <- function(x, densityNames, max.it=500, eps=1e-4) {

  # number of states
  N = length(densityNames)
  # check that all densities are in Other, Z, NB
  enum = c("Other", "Z", "NB")
  densityNames = factor(densityNames, levels=enum)
  if (any(is.na(densityNames))) {
    stop(paste("densityNames must be in", enum))
  }
  # transform to integers
  densityNames = as.numeric(densityNames) - 1
  
  # length of observation sequence
  T = length(x) 
 
  # void R_univariate_hmm(double* O, int* T, int* N, int* densityNames, double* r, double* p, double* w, int* iterationMAX, double* eps, double* post, double* A, double* proba, double *loglik)
  z = .C("R_univariate_hmm",
    as.double(x), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.integer(densityNames), # int* densityNames
    double(length=N), # double* r
    double(length=N), # double* p
    double(length=N), # double* w
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N), # double* proba
    double(length=1)) #double* loglik
  # POST will hold the results

  posterior = matrix(z[[10]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[11]], ncol=N), proba=z[[12]])
  params[["distributions"]] = cbind(r=z[[5]], p=z[[6]], w=z[[7]])
  attr(posterior, "params") = params
  attr( posterior, "loglik") = z[[13]] #TODO: check whether this works
  return(posterior)
}
