#multivariate.hmm <- function(x, Nmod, N, states=NULL, r=NULL, p=NULL, w=NULL, cor_matrix_inv=NULL, det=NULL, max.it=500, eps=1e-4) {
multivariate.hmm <- function(reads, correlationMatrixInverse, deter, params, states, max.it=500, eps=1e-4){
  #number of states
  N = length(states)
  #observation length
  T = nrow(reads)
  #number of modifications
  Nmod = ncol(reads)
  
  #get r, p, and w
  r <- c()
  p <- c()
  w <- c()
  for(imod in 1:Nmod)
  {
    for (um in 1:2)
    {
	r <- c(r, params[[names(params)[imod]]]$distributions[um,"r"])
        p <- c(p, params[[names(params)[imod]]]$distributions[um,"p"])
        w <- c(w, params[[names(params)[imod]]]$distributions[um,"w"])
    }
  }

  #void R_multivariate_hmm(double* O, int* T, int* N, int *Nmod, int* states, double* r, double* p, double* w, double* cor_matrix_inv, double* det, int* iterationMAX, double* eps, double* post, double* A, double* proba, double* loglik)
  z = .C("R_multivariate_hmm",
    as.double(as.vector(t(reads))), # double* O
    as.integer(T), # int* T
    as.integer(N), # int* N
    as.integer(Nmod), #int* Nmod
    as.integer(states), #int* states
    as.double(r), # double* r
    as.double(p), # double* p
    as.double(w), # double* w
    as.double(correlationMatrixInverse), #as.double(as.vector(correlationMatrixInverse)), #as.double(correlationMatrixInverse), #double* cor_matrix_inv
    as.double(deter), #double* det
    as.integer(max.it), #  int* iterationMAX
    as.double(eps), # double* eps
    double(length=T * N), # double* post
    double(length=N*N), # double* A
    double(length=N), # double* proba
    double(length=1))#double* loglik
  # POST will hold the results
  posterior = matrix(z[[13]], ncol=N)
  colnames(posterior) = paste("Pr(q_i=", 1:ncol(posterior)-1, "|O,model)", sep="")
  params = list(A=matrix(z[[14]], ncol=N), proba=z[[15]])
  params[["distributions"]] = cbind(r=z[[6]], p=z[[7]], w=z[[8]])
  attr(posterior, "params") = params
  attr(posterior, "loglik") = z[[16]]
#  posterior$loglik = z[[16]] ##putting the loglik into the results as well
  return(posterior)
}
