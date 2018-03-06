#include "densities.h"
#include "utility.h"


//double const pi;

ostream& operator <<(ostream &os, Density &obj) {
  obj.put(os);
  return(os);
};

void print_vector(double* x, int p) {
  cout<<"c(";
  for (int i=0; i<p-1; i++) {
    cout<<x[i]<<", ";
  }
  cout<<x[p-1]<<")";
}

void print_ivector(int* x, int p) {
  cout<<"c(";
  for (int i=0; i<p-1; i++) {
    cout<<x[i]<<", ";
  }
  cout<<x[p-1]<<")";
}

Zero::Zero(double* observations) {
  this->O = observations;
}

double Zero::density(int t) {
  return(density(this->O[t]));
}

double Zero::density(double* x) {
  if (*x == 0) {
    return(1);
  } else {
    return(0);
  }
}

void Zero::initialize(double* weight, int T) {}

void Zero::put(std::ostream& os) {
  os<<"Zero();";
}

void Zero::copy(Density* other) {}

void Zero::update(double* weight, int T) {}

double Zero::logdensity(int t) {
  return(log(density(t)));
}

double Zero::logdensity(double* x) {
  return(log(density(x)));
}

double Zero::CDF(double* x) {
  cout<<"Warning: Zero::CDF not implemented yet!\n";
  return(0);
}


double Zero::logCDF(double* x) {
  return(log(CDF(x)));
}

DensityName Zero::getType() {
  return(Z);
}

double Zero::getMean() {
  return(0);
}


LogNormal::LogNormal(double* observations, int n, double mu, double sigma) {
  this->O = observations;
  this->n = n;
  this->mu = mu;
  this->sigma = sigma;
};

double LogNormal::density(double* x) {
  return(1.0 / (x[0] * sqrt(pow(sigma, 2) * 2.0 * pi)) * exp(-pow(log(x[0]) - mu, 2) / (2.0 * pow(sigma, 2))));
}

double LogNormal::density(int t) {
  return(this->density(this->O[t]));
};

double LogNormal::logdensity(int t) {
  return(0);
}

double LogNormal::logdensity(double* x) {
  return(0);
}



double LogNormal::CDF(double* x) {
  cout<<"LogNormal::CDF not implemented yet!\n";
  return(0);
};

double LogNormal::logCDF(double* x) {
  return(log(CDF(x)));
}


// initialize the parameters with Maximum likelihood estimates
void LogNormal::initialize(double* weight, int T) {
  // ML estimate for the mean
  double mu = 0;
  for (int i=0; i<n; i++) {
    mu += log(O[i]);
  }
  mu = mu / n;
  this->mu = mu;

  // ML estimate for the sd
  double sd = 0;
  for (int i=0; i<n; i++) {
    sd += pow(log(O[i]) - mu, 2);
  }
  sd = sd / n;
  this->sigma = sd;
};


void LogNormal::put(std::ostream& os) {
  os<<"Lognormal(mu = "<<this->mu<<", sigma = "<<this->sigma<<");";
};

void LogNormal::copy(Density* other) {
  this->mu = ((LogNormal*)other)->mu;
  this->sigma = ((LogNormal*)other)->sigma;
};

void LogNormal::update(double* weight, int T) {
  int t;

  // effective size
  double n = 0;
  for (t=0; t<T; t++) {
    n += weight[t];
    // cout<<O[t]<<", "<<weight[t]<<", "<<log(O[t])<<", "<<log(O[t]) * weight[t]<<"\n";
  }

  double mu = 0;
  for (t=0; t<T; t++) {
    mu += (log(O[t]) * weight[t]);
  }
  mu = mu / n;

  double sd = 0;
  for (t=0; t<T; t++) {
    sd += (weight[t] * pow(log(O[t]) - mu, 2));
  }
  sd = sqrt(sd / n);

  this->mu = mu;
  this->sigma = sd;
};

DensityName LogNormal::getType() {
  return(Other);
}

double LogNormal::getMean() {
  return(0);
}


BivariateLogNormal::BivariateLogNormal(double* observations, int n, double* mu, double* sigma, double rho) {
  this->O = observations;
  this->n = n;
  this->mu = mu;
  this->sigma = sigma;
  this->rho = rho;

  for (int i=0; i<2; i++) {
    cout<<"mu["<<i<<"] = "<<mu[i]<<", sigma["<<i<<"] = "<<sigma[i]<<"\n";
  }

};

double BivariateLogNormal::density(int t) {
  double* x = &(this->O[2 * t]);
  return(this->density(x));
}

double BivariateLogNormal::density(double* xx) {
  /* local variables for better readibility of the density function */
  double x = xx[0];
  double mu_x = this->mu[0];
  double sigma_x = this->sigma[0];
  double y = xx[1];
  double mu_y = this->mu[1];
  double sigma_y = this->sigma[1];

  return(
	 exp(-1 / (2 * (1 - pow(rho, 2))) * 
	     (pow((log(x) - mu_x) / sigma_x, 2) - 
	      2 * rho * (log(x) - mu_x) / sigma_x * (log(y) - mu_y) / sigma_y +
	      pow((log(y) - mu_y) / sigma_y, 2))) *
	 (1 / (2 * pi * x * y * sigma_x * sigma_y * sqrt(1 - pow(rho, 2))))
	 );
};


double BivariateLogNormal::logdensity(int t) {
  return(0);
}

double BivariateLogNormal::logdensity(double* x) {
  return(0);
}



double BivariateLogNormal::CDF(double* x) {
  cout<<"BivariateLogNormal::CDF not implemented yet!\n";
  return(0);
};

double BivariateLogNormal::logCDF(double* x) {
  return(log(CDF(x)));
}



// ML estimates.. todo
void BivariateLogNormal::initialize(double* weight, int T) {
  
}

void BivariateLogNormal::put(ostream &os) {
  os<<"BivariateLognormal(mu = c("<<this->mu[0]<<", "<<this->mu[1]<<"), sigma = c("<<this->sigma[0]<<", "<<this->sigma[1]<<"));";
}

void BivariateLogNormal::copy(Density* other) {
  for (int i=0; i<2; i++) {
    this->mu[i] = ((BivariateLogNormal*)other)->mu[i];
    this->sigma[i] = ((BivariateLogNormal*)other)->sigma[i];
  }
};

void BivariateLogNormal::update(double* weight, int T) {
  
};

DensityName BivariateLogNormal::getType() {
  return(Other);
}

double BivariateLogNormal::getMean() {
  return(0);
}


BivariateNormal::BivariateNormal(double* observations, int n, double* mu, double* sigma, double rho) {
  this->O = observations;
  this->n = n;
  this->mu_x = mu[0];
  this->sigma_x = sigma[0];
  this->mu_y = mu[1];
  this->sigma_y = sigma[1];
  this->rho = rho;

};

BivariateNormal::BivariateNormal(double* observations, double mu_x, double sigma_x, double mu_y, double sigma_y, double rho) {
  this->O = observations;
  this->n = 0;
  this->mu_x = mu_x;
  this->sigma_x = sigma_x;
  this->mu_y = mu_y;
  this->sigma_y = sigma_y;
  this->rho = rho;
};


double BivariateNormal::density(int t) {
  return(density(&(this->O[2 * t])));
};

double BivariateNormal::density(double* o) {
  /* local variables for better readibility of the density function */
  double x = o[0];
  double y = o[1];

  return(
	 exp(-1 / (2 * (1 - pow(rho, 2))) * 
	     (pow((x - mu_x) / sigma_x, 2) - 
	      2 * rho * (x - mu_x) / sigma_x * (y - mu_y) / sigma_y +
	      pow((y - mu_y) / sigma_y, 2))) *
	 (1 / (2 * pi * sigma_x * sigma_y * sqrt(1 - pow(rho, 2))))
	 );
};


double BivariateNormal::logdensity(int t) {
  return(0);
}

double BivariateNormal::logdensity(double* x) {
  return(0);
}



double BivariateNormal::CDF(double* x) {
  cout<<"BivariateNormal::CDF not implemented yet!\n";
  return(0);
};

double BivariateNormal::logCDF(double* x) {
  return(log(CDF(x)));
}


// ML estimates.. todo
void BivariateNormal::initialize(double* weight, int T) {
  
}

void BivariateNormal::put(ostream &os) {
  os<<"Bivariatenormal(mu = c("<<this->mu_x<<", "<<this->mu_y<<"), sigma = c("<<this->sigma_x<<", "<<this->sigma_y<<"));";
}

void BivariateNormal::copy(Density* other) {
  this->mu_x = ((BivariateNormal*)other)->mu_x;
  this->sigma_x = ((BivariateNormal*)other)->sigma_x;
  this->mu_y = ((BivariateNormal*)other)->mu_y;
  this->sigma_y = ((BivariateNormal*)other)->sigma_y;
};

void BivariateNormal::update(double* weight, int T) {
  
};

DensityName BivariateNormal::getType() {
  return(Other);
}

double BivariateNormal::getMean() {
  return(0);
}

// linear algebra helper functions for the multivariate gaussian

// this is just a wrapper for the DGEMM function from the blas library 
void matrixMultiply(int transposeA, int transposeB, 
		    int nrowA, int ncolA, int nrowB, int ncolB, 
		    double *a, double *b, double *c) {
  /*
  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     .. Scalar Arguments ..
     DOUBLE PRECISION ALPHA,BETA
     INTEGER K,LDA,LDB,LDC,M,N
     CHARACTER TRANSA,TRANSB
     ..
     .. Array Arguments ..
     DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
     ..

  Purpose
  =======

  DGEMM  performs one of the matrix-matrix operations

     C := alpha*op( A )*op( B ) + beta*C,

  where  op( X ) is one of

     op( X ) = X   or   op( X ) = X**T,

  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

  Arguments
  ==========

  TRANSA - CHARACTER*1.
           On entry, TRANSA specifies the form of op( A ) to be used in
           the matrix multiplication as follows:

              TRANSA = 'N' or 'n',  op( A ) = A.

              TRANSA = 'T' or 't',  op( A ) = A**T.

              TRANSA = 'C' or 'c',  op( A ) = A**T.

           Unchanged on exit.

  TRANSB - CHARACTER*1.
           On entry, TRANSB specifies the form of op( B ) to be used in
           the matrix multiplication as follows:

              TRANSB = 'N' or 'n',  op( B ) = B.

              TRANSB = 'T' or 't',  op( B ) = B**T.

              TRANSB = 'C' or 'c',  op( B ) = B**T.

           Unchanged on exit.

  M      - INTEGER.
           On entry,  M  specifies  the number  of rows  of the  matrix
           op( A )  and of the  matrix  C.  M  must  be at least  zero.
           Unchanged on exit.

  N      - INTEGER.
           On entry,  N  specifies the number  of columns of the matrix
           op( B ) and the number of columns of the matrix C. N must be
           at least zero.
           Unchanged on exit.

  K      - INTEGER.
           On entry,  K  specifies  the number of columns of the matrix
           op( A ) and the number of rows of the matrix op( B ). K must
           be at least  zero.
           Unchanged on exit.

  ALPHA  - DOUBLE PRECISION.
           On entry, ALPHA specifies the scalar alpha.
           Unchanged on exit.

  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
           part of the array  A  must contain the matrix  A,  otherwise
           the leading  k by m  part of the array  A  must contain  the
           matrix A.
           Unchanged on exit.

  LDA    - INTEGER.
           On entry, LDA specifies the first dimension of A as declared
           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
           LDA must be at least  max( 1, m ), otherwise  LDA must be at
           least  max( 1, k ).
           Unchanged on exit.

  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
           part of the array  B  must contain the matrix  B,  otherwise
           the leading  n by k  part of the array  B  must contain  the
           matrix B.
           Unchanged on exit.

  LDB    - INTEGER.
           On entry, LDB specifies the first dimension of B as declared
           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
           LDB must be at least  max( 1, k ), otherwise  LDB must be at
           least  max( 1, n ).
           Unchanged on exit.

  BETA   - DOUBLE PRECISION.
           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
           supplied as zero then C need not be set on input.
           Unchanged on exit.

  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
           Before entry, the leading  m by n  part of the array  C must
           contain the matrix  C,  except when  beta  is zero, in which
           case C need not be set on entry.
           On exit, the array  C  is overwritten by the  m by n  matrix
           ( alpha*op( A )*op( B ) + beta*C ).

  LDC    - INTEGER.
           On entry, LDC specifies the first dimension of C as declared
           in  the  calling  (sub)  program.   LDC  must  be  at  least
           max( 1, m ).
           Unchanged on exit. 
  */

  int m = nrowA;
  int n = ncolB;
  int k = ncolA;
  double one = 1;
  double zero = 0;

  char transa = 'n';
  if (transposeA) {
    transa = 't';
    m = ncolA;
    k = nrowA;
  }
  char transb = 'n';
  if (transposeB) {
    transb = 't';
    n = nrowB;
    // sanity check for the dimensions
    if (k != ncolB) {
      cout<<"dimension do not match "<<k<<" != "<<ncolB<<"\n";
    }
  }
  if (k != nrowB) {
    cout<<"dimension do not match "<<k<<" != "<<nrowB<<"\n";
  }
  F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &one, a, &nrowA, b, &nrowB, &zero, c, &m);
}
  

void matrixInverse(double* x, int nrow, int ncol) {

  /* here we have to use two lapack routines 
     1) factorization using DGETRF
     2) inversion using DGETRI */
  
  /*  SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

  -- LAPACK routine (version 3.2) --
  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
     November 2006

     .. Scalar Arguments ..
     INTEGER            INFO, LDA, M, N
     ..
     .. Array Arguments ..
     INTEGER            IPIV( * )
     DOUBLE PRECISION   A( LDA, * )
     ..

  Purpose
  =======

  DGETRF computes an LU factorization of a general M-by-N matrix A
  using partial pivoting with row interchanges.

  The factorization has the form
     A = P * L * U
  where P is a permutation matrix, L is lower triangular with unit
  diagonal elements (lower trapezoidal if m > n), and U is upper
  triangular (upper trapezoidal if m < n).

  This is the right-looking Level 3 BLAS version of the algorithm.

  Arguments
  =========

  M       (input) INTEGER
          The number of rows of the matrix A.  M >= 0.

  N       (input) INTEGER
          The number of columns of the matrix A.  N >= 0.

  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).

  IPIV    (output) INTEGER array, dimension (min(M,N))
          The pivot indices; for 1 <= i <= min(M,N), row i of the
          matrix was interchanged with row IPIV(i).

  INFO    (output) INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                has been completed, but the factor U is exactly
                singular, and division by zero will occur if it is used
                to solve a system of equations.
  */

  int* ipiv = (int*)calloc(nrow, sizeof(int));
  int info = 0;
  F77_CALL(dgetrf)(&nrow, &ncol, x, &nrow, ipiv, &info);
  if (info != 0) {
    cout<<"problem with the matrix factorization\n";
  }

  /*   SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       
  -- LAPACK routine (version 3.2) --
  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
     November 2006

     .. Scalar Arguments ..
     INTEGER            INFO, LDA, LWORK, N
     ..
     .. Array Arguments ..
     INTEGER            IPIV( * )
     DOUBLE PRECISION   A( LDA, * ), WORK( * )
     ..

  Purpose
  =======

  DGETRI computes the inverse of a matrix using the LU factorization
  computed by DGETRF.

  This method inverts U and then computes inv(A) by solving the system
  inv(A)*L = inv(U) for inv(A).

  Arguments
  =========

  N       (input) INTEGER
          The order of the matrix A.  N >= 0.

  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the factors L and U from the factorization
          A = P*L*U as computed by DGETRF.
          On exit, if INFO = 0, the inverse of the original matrix A.

  LDA     (input) INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).

  IPIV    (input) INTEGER array, dimension (N)
          The pivot indices from DGETRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i).

  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.

  LWORK   (input) INTEGER
          The dimension of the array WORK.  LWORK >= max(1,N).
          For optimal performance LWORK >= N*NB, where NB is
          the optimal blocksize returned by ILAENV.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.

  INFO    (output) INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                singular and its inverse could not be computed.

  */
  
  int* order = (int*)calloc(2, sizeof(int));
  order[0] = nrow;
  order[1] = ncol;

  double qwork;
  int lwork = -1;
  // query optimal size of the work array and allocate the memory
  F77_CALL(dgetri)(order, x, &nrow, ipiv, &qwork, &lwork, &info);
  lwork = (int)qwork;
  double* work = (double*)calloc(lwork, sizeof(double));

  // run the subroutine
  F77_CALL(dgetri)(order, x, &nrow, ipiv, work, &lwork, &info);
  if (info != 0) {
    cout<<"problem with the matrix inverse\n";
  }
  
  // free the memory
  free(work);
  free(order);
}


double determinant(double* x, int nrow, int ncol) {
  // use the R function (maybe this is easier than using lapack directly)
  SEXP A = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
  memcpy((void*)REAL(A), (void*)x, sizeof(double) * nrow * ncol);
  SEXP LOG = PROTECT(Rf_allocVector(LGLSXP, 1));
  LOGICAL(LOG)[0] = 0;

  SEXP e;
  SEXP f;
  // s will be an expression for the function setseed
  PROTECT(e = f = Rf_allocList(3));    // this is some kind of linked list
  SET_TYPEOF(e, LANGSXP);
  SETCAR(f, Rf_install("determinant"));    // CAR is the current value
  f = CDR(f);                       // CDR is a pointer to the next value
  SETCAR(f, A); 
  f = CDR(f);
  SETCAR(f, LOG);
  SET_TAG(f, Rf_install("logarithm"));
  SEXP res = Rf_eval(e, R_GlobalEnv);
  double det = REAL(VECTOR_ELT(res, 0))[0];
  UNPROTECT(3);
  return(det);
}


MultivariateNormal::MultivariateNormal(double* observations, double* mu, double* sigma, int p) {
  this->O = observations;
  this->mu = mu;
  this->p = p;
  this->sigma_inv = (double*) calloc(p * p, sizeof(double));
  setSigma(sigma);
}

MultivariateNormal::~MultivariateNormal() {
  free(this->sigma_inv);
}

double MultivariateNormal::density(int t) {
  double* x = &(this->O[t * this->p]);
  return(density(x));
}

double MultivariateNormal::density(double* x) { 
  /* we do matrix multiplication using the blas library check the comments in
     the matrixMultiply function which is a wrapper for DGEMM
  */
  
  // compute (x - mu)^t Sigma^-1 (x - mu)
  double* z = (double*) calloc(p, sizeof(double));
  for (int i=0; i<p; i++) {
    z[i] = x[i] - mu[i];
  }

  // compute z^t Sigma^-1 z
  // left part
  double* left = (double*)calloc(p, sizeof(double));
  matrixMultiply(1, 0, p, 1, p, p, z, sigma_inv, left);
  
  // right part
  double exponent = 0;
  matrixMultiply(0, 0, 1, p, p, 1, left, z, &exponent);
  free(left);
  exponent = -0.5 * exponent;

  double d = 1 / (pow(2 * pi, ((double)p) / 2.0) * sqrt(sigma_det)) * exp(exponent);
  return(d);
}

double MultivariateNormal::logdensity(int t) {
  return(0);
}

double MultivariateNormal::logdensity(double* x) {
  return(0);
}



double MultivariateNormal::CDF(double* x) {
  cout<<"MultivariateNormal::CDF not implemented yet!\n";
  return(0);
};

double MultivariateNormal::logCDF(double* x) {
  return(log(CDF(x)));
}


void MultivariateNormal::initialize(double* weight, int T) {

}

void MultivariateNormal::put(std::ostream& os) {
  os<<"MultivariateNormal(mu=c(";
  for (int i=0; i<p - 1; i++) {
    os<<mu[i]<<", ";
  }
  os<<mu[p - 1]<<"), sigma=c(";
  for (int i=0; i<p * p - 1; i++) {
    os<<this->sigma[i]<<", ";
  }
  os<<this->sigma[p * p - 1]<<"));\n";
}

void MultivariateNormal::copy(Density* other) {
  
}

void MultivariateNormal::update(double* weight, int T) {

}

// take care to compute the inverse and determinant so density computation
// can use this without recomputing
void MultivariateNormal::setSigma(double* sigma) {
  this->sigma = sigma;
  // compute the inverse of sigma
  memcpy((void*)sigma_inv, (void*)(this->sigma), sizeof(double) * p * p);
  matrixInverse(sigma_inv, p, p);
  // compute the determinant of sigma
  sigma_det = determinant(this->sigma, p, p);
}

DensityName MultivariateNormal::getType() {
  return(Other);
}

double MultivariateNormal::getMean() {
  return(0);
}


/* for testing of the mv normal density */
extern "C" {
  void R_dmvn(double* O, int* N, int* p, double* mu, double* sigma, double* density) {
    MultivariateNormal* mvn = new MultivariateNormal(O, mu, sigma, *p);

    cout<<*mvn<<"\n";

    for (int i=0; i<*N; i++) {
      density[i] = mvn->density(i);
    }
  }
}


Zinba::Zinba(double* observations, int n, double size, double mu, double beta) {
  this->O = observations;
  this->n = n;
  this->size = size;
  this->mu = mu;
  this->beta = beta;
  this->max_x = 0;
  this->logcdf = new vector<double>();
  // initialize
  double zero = 0;
  logcdf->push_back(logdensity(&zero));
}

Zinba::Zinba(double size, double mu, double beta) {
  this->O = NULL;
  this->n = 0;
  this->size = size;
  this->mu = mu;
  this->beta = beta;
  this->max_x = 0;
  this->logcdf = new vector<double>();
  // initialize
  double zero = 0;
  logcdf->push_back(logdensity(&zero));
}

double Zinba::density(double x, double size, double mu, double beta) {
  double prob = size / (size + mu);
  // last param is log
  double density = (1.0 - beta) * dnbinom(x, size, prob, 0); 
  if (x == 0) {
    density += beta;
  }
  return(density);
}

double Zinba::density(double* x) {
  return(density(*x, size, mu, beta));
}

double Zinba::density(int t) {
  return(density(O[t], size, mu, beta));
}

double Zinba::logdensity(int t) {
  return(logdensity(O[t]));
}

double Zinba::logdensity(double* x) {
  double prob = size / (size + mu);
  // last param is log
  double density = log(1.0 - beta) + dnbinom(*x, size, prob, 1); 
  if (x == 0) {
    density = logspace_add(log(beta), density);
  }
  return(density);
}


double Zinba::logCDF(double* x) {
  int k = (int) *x;
  if (k < 0) {
    return(1);
  }
  // check if we have already computed the cdf this far
  if (k <= max_x) {
    double value = this->logcdf->at(k);
    return(value);
  } else {
    // if not fill up the missing parts between max_x and x
    double current = 0;
    // switch the parameterization of the negative binomial to use pbeta
    double p = size / (size + mu);
    for (int i=max_x+1; i<=k; i++) {
      // here we actually compute the cdf in log space
      // last 2 args of pbeta are lower tail and log
      current = logspace_add(log(beta), log(1 - beta) + pbeta(p, size, (double)i + 1, 1, 1));
      logcdf->push_back((double) current);
      // cout<<i<<", "<<d<<", "<<log(d)<<", "<<current<<", "<<log(current)<<"\n";
    }
    max_x = k;
    return(current);
  }

}

double Zinba::CDF(double* x) {
  return(exp(logCDF(x)));
}


void Zinba::initialize(double* weight, int T) {
  // beta is estimated just as the probortion of zeros
  int nzero = 1; // make sure it is not zero
  for (int i=0; i<n; i++) {
    if (O[i] == 0) {
      nzero++;
    }
  }
  beta = (double)nzero / (double)(n + 1);
  // this is taken from the fitdistr function in MASS
  mu = 0;
  for (int i=0; i<n; i++) {
    mu += O[i];
  }
  mu = mu / n;
  double v = 0;
  for (int i=0; i<n; i++) {
    v += pow(O[i] - mu, 2);
  }
  v = v / n;
  if (v > mu) {
    size = pow(mu, 2) / (v - mu);
  } else {
    size = 100;
  }
}

void Zinba::initialize(double* weight) {
  // beta is estimated just as the probortion of zeros
  double sum = 0;
  double nzero = 0;
  for (int i=0; i<n; i++) {
    sum += weight[i];
    if (O[i] == 0) {
      nzero += weight[i];
    }
  }
  beta = (double)nzero / sum;
  // this is taken from the fitdistr function in MASS
  mu = 0;
  for (int i=0; i<n; i++) {
    mu += weight[i] * O[i];
  }
  mu = mu / sum;
  double v = 0;
  for (int i=0; i<n; i++) {
    v += weight[i] * pow(O[i] - mu, 2);
  }
  v = v / sum;
  if (v > mu) {
    size = pow(mu, 2) / (v - mu);
  } else {
    size = 100;
  }
}



void Zinba::put(ostream &os) {
  os<<"Zinba(mu = "<<this->mu<<", size = "<<this->size<<", beta = "<<this->beta<<");";
}

void Zinba::copy(Density* other) {
  this->mu = ((Zinba*)other)->mu;
  this->size = ((Zinba*)other)->size;
  this->beta = ((Zinba*)other)->beta;
}

double Zinba::logLikelihood() {
  double loglik = 0;
  for (int i=0; i<n; i++) {
    loglik += log(this->weight[i]) + this->logdensity(i);
  }
  return(loglik);
}

double Zinba_objective(int n, double *par, void *ex) {
  // this is the objective function to maximize the expected log likelihood
  // ex is a pointer to the original Zinba object 

  // check if parameters are valid
  if (par[0] < 0 || par[1] < 0 || par[2] < 0 || par[2] > 1) {
    return(INFINITY);
  }

  Zinba* obj = (Zinba*)ex;
  obj->mu = par[0];
  obj->size = par[1];
  obj->beta = par[2];

  return(-(obj->logLikelihood())); // optim is minimizing..
}

void Zinba::update(double* weight, int T) {
  // assign the weights to the temp storage (actually just a pointer..)
  this->weight = weight;

  // starting parameters are the current parameters
  double* xin = (double*) calloc(3, sizeof(double));
  xin[0] = mu;
  xin[1] = size;
  xin[2] = beta;

  // alloc space for the final params
  double* x = (double*) calloc(3, sizeof(double));

  // compute the current loglikelihood to check if optimization improves
  double oldLogLik = logLikelihood();
  Zinba* old = new Zinba(size, mu, beta);

  // the new log likelihood is the final result of the optimization
  double newLogLik;

  // status can be 0: OK, 1: maxit reached, 10: degenerate simplex
  int fail;

  // tolerance (abs is -Inf, and relative is machine prec. from R source optim)
  double intol = sqrt(DBL_EPSILON);

  // ex is the extra argument passed to the objective fn (this Zinba object)
  void* ex = (void*)this;

  // counter for function evaluations
  int fncount = 0;

  int maxit = 100;

  // defaults for alpha=1, beta=0.5 and gamma=2 from the R manual
  // use the Nelder Mead algorithm from R_ext/Applic.h
  nmmin(3, xin, x, &newLogLik, Zinba_objective, &fail, -INFINITY, intol, ex,
        1.0, 0.5, 2.0, 0, &fncount, maxit);
  newLogLik = -newLogLik;

  // if the new loglikelihood is smaller than the previous we stick to the
  // old parameters
  if (newLogLik <= oldLogLik) {
    cout<<"optimization failed, reverting back to old parameters\n";
    cout<<"old: "<<*old<<"\n";
    this->copy(old);
  }

  cout<<"previous loglik "<<oldLogLik<<" new loglik "<<newLogLik<<" delta "<<oldLogLik - newLogLik<<" nelder mead status: "<<fail<<"\n";

  free(xin);
  free(x);
};


DensityName Zinba::getType() {
  return(Other);
}

double Zinba::getMean() {
  return(0);
}


// testing this in R
extern "C" {
  void R_fitzinba(double* x, double* w, int* n, double* mu, double* size, double* beta) {
    Zinba* zinba = new Zinba(x, *n, *size, *mu, *beta);
    zinba->initialize(w, *n);
    cout<<"Initialized: "<<*zinba<<"\n";
    zinba->update(w, *n);
    cout<<"Fitted: "<<*zinba<<"\n";
    *mu = zinba->mu;
    *size = zinba->size;
    *beta = zinba->beta;
  }
}


Negbinom::Negbinom(double* observations, int n, double size, double mu) {
  this->O = observations;
  this->n = n;
  this->size = size;
  this->mu = mu;
  this->max_x = 0;
  this->logcdf = new vector<double>();
  // initialize
  logcdf->push_back(log(density(0, size, mu)));
}

Negbinom::Negbinom(double size, double mu) {
  this->O = NULL;
  this->n = 0;
  this->size = size;
  this->mu = mu;
  this->max_x = 0;
  this->logcdf = new vector<double>();
  // initialize
  logcdf->push_back(log(density(0, size, mu)));
}

double Negbinom::density(double x, double size, double mu) {
  double prob = size / (size + mu);
  // last param is log
  double density = dnbinom(x, size, prob, 0); 
  return(density);
}

double Negbinom::logdensity(double* x) {
  double prob = size / (size + mu);
  // last param is log
  double density = dnbinom(*x, size, prob, 1); 
  return(density);
}

double Negbinom::logdensity(int t) {
  return(logdensity(&(O[t])));
}

double Negbinom::density(int t) {
  return(density(O[t], size, mu));
}

double Negbinom::density(double* x) {
  return(density(*x, size, mu));
}


double Negbinom::logCDF(double* x) {
  int k = (int) *x;
  if (k < 0) {
    return(1);
  }
  // check if we have already computed the cdf this far
  if (k <= max_x) {
    double value = this->logcdf->at(k);
    return(value);
  } else {
    // if not fill up the missing parts between max_x and x
    double current = 0;
    // switch the parameterization of the negative binomial to use pbeta
    double p = size / (size + mu);
    for (int i=max_x+1; i<=k; i++) {
      // here we actually compute the cdf in log space
      // last 2 args of pbeta are lower tail and log
      current = pbeta(p, size, (double)i + 1, 1, 1);
      logcdf->push_back((double) current);
      // cout<<i<<", "<<d<<", "<<log(d)<<", "<<current<<", "<<log(current)<<"\n";
    }
    max_x = k;
    return(current);
  }

}

double Negbinom::CDF(double* x) {
  return(exp(logCDF(x)));
}


void Negbinom::initialize() {
  // this is taken from the fitdistr function in MASS
  mu = 0;
  for (int i=0; i<n; i++) {
    mu += O[i];
  }
  mu = mu / n;
  double v = 0;
  for (int i=0; i<n; i++) {
    v += pow(O[i] - mu, 2);
  }
  v = v / n;
  if (v > mu) {
    size = pow(mu, 2) / (v - mu);
  } else {
    size = 100;
  }
}

void Negbinom::initialize(double* weight, int T) {
  // this is taken from the fitdistr function in MASS
  double sum = 0;
  mu = 0;
  for (int i=0; i<n; i++) {
    sum += weight[i];
    mu += weight[i] * O[i];
  }
  mu = mu / sum;
  double v = 0;
  for (int i=0; i<n; i++) {
    v += weight[i] * pow(O[i] - mu, 2);
  }
  v = v / sum;
  if (v > mu) {
    size = pow(mu, 2) / (v - mu);
  } else {
    size = 100;
  }
}

DensityName Negbinom::getType() {
  return(NB);
}

double Negbinom::getMean() {
  return(mu);
}

void Negbinom::put(ostream &os) {
  os<<"Negbinom(mu = "<<this->mu<<", size = "<<this->size<<");";
}

void Negbinom::copy(Density* other) {
  this->mu = ((Negbinom*)other)->mu;
  this->size = ((Negbinom*)other)->size;
}


double Negbinom::logLikelihood() {
  double logLik = 0;
  for (int i=0; i<n; i++) {
    logLik += weight[i] * log(this->density(i)); // this will be adapted
  }
  return(logLik);
}

double Negbinom::logLikelihoodPartialDerivative(double size, void* extra) {
  Negbinom* self = (Negbinom*) extra;
  // partial derivative of the loglikelihood function with respect to size
  // substituting already the ML soluation for p / mu (see update function)
  double l = 0;
  double sum_w = 0;
  double sum_wy = 0;
  for (int i=0; i<self->n; i++) {
    sum_w += self->weight[i];
    sum_wy += self->weight[i] * self->O[i];
  }
  double p = size * sum_w / (sum_wy + size * sum_w);

  for (int i=0; i<self->n; i++) {
    l += self->weight[i] * (digamma(self->O[i] + size) - digamma(size) + log(p));
  }
  // cout<<"partial derivative f("<<size<<") = "<<l<<"\n";
  return(l);
}


// code from maria to compare to
void Negbinom::maria_update(double* weight, int T) {
    int ti, jN, k, kmax;
    double eps = 1e-4;
    double numerator, denominator, rhere,dr,Fr,dFrdr,DigammaR,DigammaRplusDR,DigammaRplusX,DigammaRplusDRplusX;

    double r = this->size;

    ///////////////// UPDATE p
    numerator=denominator=0.0;
    for(ti=0;ti<T;ti++){
        numerator+=weight[ti]*r;
        denominator+=weight[ti]*(r+this->O[ti]);
    };

    double phere = numerator/denominator;//////CAUTION! THE UPDATE OF r IS NOW DONE WITH THIS UPDATED p
    ///////////////// UPDATE r
    ///////////////   NEWTON METHOD
    ///////////////////////////////////////////////////////////////////////
	rhere=r;
	dr=0.00001;
	kmax=20;
    for(k=1;k<kmax;k++){
        Fr=dFrdr=0.0;
        DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
        DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
        for(int ti=0;ti<T;ti++){
	  DigammaRplusX = digamma(rhere+this->O[ti]); //boost::math::digamma<>(rhere+this->O[ti]);
	  DigammaRplusDRplusX = digamma(rhere+dr+this->O[ti]); // boost::math::digamma<>(rhere+dr+this->O[ti]);
            if(this->O[ti]==0){
                Fr+=weight[ti]*log(phere);
                //dFrdr+=0;
            };
            if(this->O[ti]!=0){
                Fr+=weight[ti]*(log(phere)-DigammaR+DigammaRplusX);
                dFrdr+=weight[ti]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
            };
        };
        if(fabs(Fr)<eps) break;
        if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
        if(Fr/dFrdr>rhere) rhere=rhere/2.0;
        if(k%50==0){
            cout<<"NewtonRaphson: k= "<<k<<", Fr= "<<Fr<<", dr="<<-Fr/dFrdr<<", r="<<rhere<<endl;
        };
    };
    double mu = rhere * (1 - phere) / phere;
    cout<<"Maria update: r= "<<rhere<<" , p= "<<phere<<" mu = "<<mu<<endl;
}

void Negbinom::update(double* weight, int T) {
  // keep the old parameters
  // Negbinom* old = new Negbinom(size, mu);

  this->weight = weight;
  this->n = T;

  double tol = 1e-4;
  // find the root of the partial derivative numerically
  size = zeroin(DBL_EPSILON, 1000, &(Negbinom::logLikelihoodPartialDerivative), tol, this);

  // update mu (use p)
  // we obtained this by setting the partial derivative of the expected 
  // loglikelihood to zero
  double sum_w = 0;
  double sum_wy = 0;
  for (int i=0; i<T; i++) {
    sum_w += weight[i];
    sum_wy += weight[i] * O[i];
  }
  double p = size * sum_w / (sum_wy + size * sum_w);
  mu = size * (1 - p) / p;
  cout<<"Negbinom::update: p = "<<p<<" mu = "<<mu<<"\n";

};


// testing this in R
extern "C" {
  void R_fitnegbinom(double* x, double* w, int* n, double* mu, double* size) {
    Negbinom* nb = new Negbinom(x, *n, *size, *mu);

    nb->update(w, *n);
    cout<<"Fitted: "<<*nb<<"\n";
    *mu = nb->mu;
    *size = nb->size;

    // for comparison run marias update function
    // nb->maria_update(w, *n);
  }
}


BivariateZinba::BivariateZinba(double* observations, int n, double marginal_size, double marginal_mu, double marginal_beta, double* size_coef, double* mu_coef, double beta0) {
  this->O = observations;
  this->n = n;
  this->marginal_size = marginal_size;
  this->marginal_mu = marginal_mu;
  this->marginal_beta = marginal_beta;
  this->size_coef = size_coef;
  this->mu_coef = mu_coef;
  this->beta0 = beta0;
}





double BivariateZinba::density(int t) {
  double* x = &(this->O[2* t]);
  return(density(x));
}
  
double BivariateZinba::density(double* o) {
  double x = o[0]; // this follows the conditional distr
  double y = o[1]; // this follows the marginal distr
  
  double size = size_coef[0] + size_coef[1] * y;
  double mu = mu_coef[0] + mu_coef[1] * y;
  
  // the beta decays exponentially
  double beta = beta0 * exp(-y);

  double p_x = Zinba::density(x, size, mu, beta);
  double p_y = Zinba::density(y, marginal_size, marginal_mu, marginal_beta);
  return(p_x * p_y);
}

double BivariateZinba::logdensity(int t) {
  return(0);
}

double BivariateZinba::logdensity(double* x) {
  return(0);
}




double BivariateZinba::CDF(double* x) {
  cout<<"BivariateZinba::CDF not implemented yet!\n";
  return(0);
};

double BivariateZinba::logCDF(double* x) {
  return(log(CDF(x)));
}


void BivariateZinba::initialize(double* weight, int T) {
  
}

void BivariateZinba::put(std::ostream& os) {
  os<<"BivariateZinba(marginal_mu = "<<this->marginal_mu<<", marginal_size = "<<this->marginal_size<<", beta0 = "<<this->beta0<<", mu_coef = c("<<this->mu_coef[0]<<", "<<this->mu_coef[1]<<"), size_coef = c("<<this->size_coef[0]<<", "<<this->size_coef[1]<<"));";
}

void BivariateZinba::copy(Density* other) {
  
}

void BivariateZinba::update(double* weight, int T) {
  
};

DensityName BivariateZinba::getType() {
  return(Other);
}

double BivariateZinba::getMean() {
  return(0);
}

ZinbaCopula::ZinbaCopula(double* O, double size_x, double mu_x, double beta_x, double size_y, double mu_y, double beta_y, double sigma_x, double sigma_y, 
double rho) {
  this->O = O;
  this->size_x = size_x;
  this->mu_x = mu_x;
  this->beta_x = beta_x;
  this->size_y = size_y;
  this->mu_y = mu_y;
  this->beta_y = beta_y;
  this->sigma_x = sigma_x;
  this->sigma_y = sigma_y;
  this->rho = rho;
 
  this->copula_mu_x = 0;
  this->copula_mu_y = 0;
 
  /* density evaluation is expensive, so we store all values that were 
     previously computed in a sparse matrix */
  
//   density evaluation is expensive, so we store all values that were previously computed in a sparse matrix 
  
  this->zinba_x = new Zinba(size_x, mu_x, beta_x);
  this->zinba_y = new Zinba(size_y, mu_y, beta_y);

  this->bvn = new BivariateNormal(O, 0, sigma_x, 0, sigma_y, rho);
}


ZinbaCopula::ZinbaCopula(double* O, double size_x, double mu_x, double beta_x, double size_y, double mu_y, double beta_y, double sigma_x, double sigma_y, double rho, double copula_mu_x, double copula_mu_y) {
  this->O = O;
  this->size_x = size_x;
  this->mu_x = mu_x;
  this->beta_x = beta_x;
  this->size_y = size_y;
  this->mu_y = mu_y;
  this->beta_y = beta_y;
  this->sigma_x = sigma_x;
  this->sigma_y = sigma_y;
  this->rho = rho;
 
  this->copula_mu_x = copula_mu_x;
  this->copula_mu_y = copula_mu_y;
 
  /* density evaluation is expensive, so we store all values that were 
     previously computed in a sparse matrix */
  
//   density evaluation is expensive, so we store all values that were previously computed in a sparse matrix 
  
  this->zinba_x = new Zinba(size_x, mu_x, beta_x);
  this->zinba_y = new Zinba(size_y, mu_y, beta_y);

  this->bvn = new BivariateNormal(O, 0, sigma_x, 0, sigma_y, rho);
}



// for the density we need to interface some fortran code which needs access
// to R's random generator, so we make it available here:
extern "C" {
  void F77_SUB(rndstart)(void) { GetRNGstate(); }
  void F77_SUB(rndend)(void) { PutRNGstate(); }
  double F77_SUB(unifrnd)(void) { return R::unif_rand(); }
}

double ZinbaCopula::density(int t) {
  double* x = &(this->O[t * 2]);
  return(density(x));
}

double ZinbaCopula::density(double* o) { 
  double x = o[0];
  double y = o[1];


  // we also want to get the same results when calling with a specific pair x,y
  // so we set R's random seed

  //  SEXP e;
  //  SEXP f;
  //  // s will be an expression for the function setseed
  //  PROTECT(e = f = allocList(2));    // this is some kind of linked list
  //  SET_TYPEOF(e, LANGSXP);
  //  SETCAR(f, install("set.seed"));    // CAR is the current value
  //  f = CDR(f);                       // CDR is a pointer to the next value
  //  SETCAR(f, ScalarReal(0)); 
  //  SET_TAG(f, install("seed"));
  //  eval(e, R_GlobalEnv);
  //  UNPROTECT(1);
  //

#ifdef DEBUG
  int print_condition = (t >= 2352 && t <= 2357 || t < 5);
  if (print_condition) {
    cout<<"ZinbaCopula::density("<<t<<") x = "<<x<<" y = "<<y<<"\n";
  }
#endif

     //check if the density for this pair has been computed already
     //precomputed values are stored in a hashmap
//  map<int, map<int, double>*>::iterator it_x;
//  map<int, double>::iterator it_y;
//  it_x = precomputed->find((int)x);
//  // if not found the iterator is equal to end
//  if (it_x != precomputed->end()) { 
//    // the hashmap values are stored in the field second (keys in first)
//    it_y = it_x->second->find((int)y);
//    if (it_y != it_x->second->end()) {
//      return(it_y->second);
//    }
//  }
//
     //compute for a new pair 

  // upper limits for the integration 
  // transform to p-values (uniform on [0, 1]) 
  double log_p_x = zinba_x->logCDF(&x);
  double log_p_y = zinba_y->logCDF(&y);

  // transform to z-scores and store them in a double pointer 
  double* upper = (double*)calloc(2, sizeof(double));
  // the last two arguments to qnorm are lower tail / logp
  upper[0] = qnorm(log_p_x, 0, 1, 1, 1) / sigma_x;
  upper[1] = qnorm(log_p_y, 0, 1, 1, 1) / sigma_y;

  // we also need the lower limits for the integration 
  double x_minus_one = x - 1;
  double y_minus_one = y - 1;
  double log_lp_x = zinba_x->logCDF(&x_minus_one);
  double log_lp_y = zinba_y->logCDF(&y_minus_one);
  double* lower = (double*)calloc(2, sizeof(double));
  lower[0] = qnorm(log_lp_x, 0, 1, 1, 1) / sigma_x;
  lower[1] = qnorm(log_lp_y, 0, 1, 1, 1) / sigma_y;

  //   cout<<"p_x = "<<exp(log_p_x)<<" p_y = "<<exp(log_p_y)<<"\n";
  //   cout<<"log(p_x) = "<<log_p_x<<" log(p_y) = "<<log_p_y<<"\n";
  //   cout<<"lp_x = "<<exp(log_lp_x)<<" lp_y = "<<exp(log_lp_y)<<"\n";
  //   cout<<"log(lp_x) = "<<log_lp_x<<" log(lp_y) = "<<log_lp_y<<"\n"; 


  /* evaluate the multivariate gaussian CDF using the following fortran routine
     (it is copied from the mvtnorm package: mvt.f)

      SUBROUTINE MVTDST( N, NU, LOWER, UPPER, INFIN, CORREL, DELTA, 
                         MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM )       

     A subroutine for computing non-central multivariate t probabilities.
     This subroutine uses an algorithm (QRSVN) described in the paper
     "Comparison of Methods for the Computation of Multivariate 
         t-Probabilities", by Alan Genz and Frank Bretz
         J. Comp. Graph. Stat. 11 (2002), pp. 950-971.

          Alan Genz 
          Department of Mathematics
          Washington State University 
          Pullman, WA 99164-3113
          Email : AlanGenz@wsu.edu

       Original source available from
       http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f

       This is version 7/10 with better support for 100 < dimension < 1000

     Parameters

     N      INTEGER, the number of variables.    
     NU     INTEGER, the number of degrees of freedom.
            If NU < 1, then an MVN probability is computed.
     LOWER  DOUBLE PRECISION, array of lower integration limits.
     UPPER  DOUBLE PRECISION, array of upper integration limits.
     INFIN  INTEGER, array of integration limits flags:
             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
     CORREL DOUBLE PRECISION, array of correlation coefficients; 
            the correlation coefficient in row I column J of the 
            correlation matrixshould be stored in 
               CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
            The correlation matrix must be positive semi-definite.
     DELTA  DOUBLE PRECISION, array of non-centrality parameters.
     MAXPTS INTEGER, maximum number of function values allowed. This 
            parameter can be used to limit the time. A sensible 
            strategy is to start with MAXPTS = 1000*N, and then
            increase MAXPTS if ERROR is too large.
     ABSEPS DOUBLE PRECISION absolute error tolerance.
     RELEPS DOUBLE PRECISION relative error tolerance.
     ERROR  DOUBLE PRECISION estimated absolute error, 
            with 99% confidence level.
     VALUE  DOUBLE PRECISION estimated value for the integral
     INFORM INTEGER, termination status parameter:
            if INFORM = 0, normal completion with ERROR < EPS;
            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
                           function vaules used; increase MAXPTS to 
                           decrease ERROR;
            if INFORM = 2, N > 1000 or N < 1.
            if INFORM = 3, correlation matrix not positive semi-definite.
  */

  int N = 2;
  int df = 0;

  // set the infinity indicators 
  int* inf = (int*)calloc(2, sizeof(int));
  if (x == 0) {
    inf[0] = 0;
  } else {
    inf[0] = 2;
  }
  // we need to be careful if p_x = 1 then the upper limit is inf 
  if (log_p_x == 0) {
    inf[0] = 1;
  }
  if (y == 0) {
    inf[1] = 0;
  } else {
    inf[1] = 2;
  }
  // we need to be careful if p_y = 1 then the upper limit is inf 
  if (log_p_y == 0) {
    inf[1] = 1;
  }


  double cor = this->rho;
  // mean vector for the copula is (0, 0) 
  double* mu = (double*)calloc(2, sizeof(double));

  // these are the defaults from the mvtnorm R package 
  int maxfuncalls = 25000;
  double releps = 0;
  double abseps = 0.001;
  double error = 0;
  int info = 0;
  double value = 0;

  //  cout<<"call to mvtdst with upper = c("<<upper[0]<<", "<<upper[1]<<") ";
  //   cout<<"lower = c("<<lower[0]<<", "<<lower[1]<<")\n";
  //   cout<<"inf = c("<<inf[0]<<", "<<inf[1]<<")\n"; 
  

  F77_CALL(mvtdst)(&N, &df, lower, upper, inf, &cor, mu, &maxfuncalls, &abseps, &releps, &error, &value, &info);
  
#ifdef DEBUG
  if (print_condition) {
    cout<<" value: "<<value<<" info: "<<info<<" error: "<<error;
  }


  if (ISNAN(value)) {
    cout<<"NAN value after mvtdst: "<<value<<" info: "<<info<<" error: "<<error<<"\n";
    cout<<"call to mvtdst with\nx = "<<x<<" y = "<<y<<"\n";
    cout<<"p_x = "<<exp(log_p_x)<<" p_y = "<<exp(log_p_y)<<"\n";
    cout<<"log_p_x = "<<log_p_x<<" log_p_y = "<<log_p_y<<"\n";
    cout<<"upper = c("<<upper[0]<<", "<<upper[1]<<")\n";
    cout<<"lower = c("<<lower[0]<<", "<<lower[1]<<")\n";
    cout<<"inf = c("<<inf[0]<<", "<<inf[1]<<")\n";
  }
#endif
  
  // there are problems with zero density values
  // if the integral gets very small, we use a rieman approximation
  // the area of the integration interval is one, so we use the normal
  // density at the upper and lower limits to find a upper bound
  if (value == 0 || ISNAN(value)) {
    double* density_in_corners = (double*)calloc(4, sizeof(double));
    density_in_corners[0] = bvn->density(upper);
    density_in_corners[1] = bvn->density(lower);
    double* corner = (double*)calloc(2, sizeof(double));
    // mix upper and lower limits to get the other two corners
    corner[0] = upper[0];
    corner[1] = lower[1];
    density_in_corners[2] = bvn->density(corner);
    corner[0] = lower[0];
    corner[1] = upper[1];
    density_in_corners[3] = bvn->density(corner);
    value = density_in_corners[0];
    for (int i=0; i<4; i++) {
      // cout<<"density(corner["<<i<<"]) = "<<density_in_corners[i]<<"\n";
      if (value < density_in_corners[i]) {
	value = density_in_corners[i];
      }
    }
    free(density_in_corners);
    free(corner);
  }

  // most likely one of the z values is Inf
  if (ISNAN(value)) {
    value = 0;
  }


  // we also observed some numerical instabilities where the value gets
  // negative even though integration limits should yield positive results
  if (value < 0) {
    // not sure what to do, just make it positive or approximate?
    value = fabs(value);
  }

  // also make sure the the upper and lower limits are actually different
  if (log_p_x == log_lp_x || log_p_y == log_lp_y) {
    value = 0;
  }

#ifdef DEBUG
  if (print_condition) {
    cout<<" final value: "<<value<<"\n";
  }
#endif

  // free the auxillary variables 
  free(lower);
  free(upper);
  free(inf);
  free(mu);
  
  // store the results for later lookup 
//  if (it_x == precomputed->end()) { 
//    map<int, double>* ymap = new map<int, double>;
//    precomputed->insert(pair<int, map<int, double>*>((int)x, ymap));
//    it_x = precomputed->find((int)x);
//  }
//  it_x->second->insert(pair<int, double>((int)y, value));
//
  return(value);
}

double ZinbaCopula::logdensity(int t) {
  return(0);
}

double ZinbaCopula::logdensity(double* x) {
  return(0);
}


  
double ZinbaCopula::CDF(double* x) {
  cout<<"ZinbaCopula::CDF not implemented yet!\n";
  return(0);
}

double ZinbaCopula::logCDF(double* x) {
  return(log(CDF(x)));
}


void ZinbaCopula::initialize(double* weight, int T) {
  
}

void ZinbaCopula::put(std::ostream& os) {
  os<<"ZinbaCopula(size_x="<<size_x<<", mu_x="<<mu_x<<", beta_x="<<beta_x<<", ";
  os<<"size_y="<<size_y<<", mu_y="<<mu_y<<", beta_y="<<beta_y<<", ";
  os<<"sigma_x="<<sigma_x<<", sigma_y="<<sigma_y<<", rho="<<rho<<");";
}

void ZinbaCopula::copy(Density* other) {
  
}

void ZinbaCopula::update(double* weight, int T) {
  
};

DensityName ZinbaCopula::getType() {
  return(Other);
}

double ZinbaCopula::getMean() {
  return(0);
}


/* for testing of the copula density */
extern "C" {
  void R_dzinbacopula(double* O, int* N, double* size_x, double* mu_x, double* beta_x, double* size_y, double* mu_y, double* beta_y, double* sigma_x, double* sigma_y, double* rho, double* copula_mu_x, double* copula_mu_y, double* density) {
    ZinbaCopula* zc = new ZinbaCopula(O, *size_x, *mu_x, *beta_x, *size_y, *mu_y, *beta_y, *sigma_x, *sigma_y, *rho, *copula_mu_x, *copula_mu_y);

    // cout<<*zc<<"\n";

    for (int i=0; i<*N; i++) {
      density[i] = zc->density(i);
    }
  }
}


MVZinbaCopula::MVZinbaCopula(double* O, double* size, double* mu, double* beta, double* sigma, int p) {
  this->O = O;
  this->size = size;
  this->mu = mu;
  this->beta = beta;
  this->p = p;

  // alloc memory for the upper triangle of the correlation matrix
  this->rho = (double*)calloc(this->p * (this->p - 1) / 2, sizeof(double));
  setSigma(sigma);

  // this->zinba = new vector<Zinba*>;
  for (int i=0; i<p; i++) {
    zinba.push_back(new Zinba(size[i], mu[i], beta[i]));
  }
  this->mvn_mu = (double*)calloc(p, sizeof(double));
  this->mvn = new MultivariateNormal(O, mvn_mu, sigma, p);
}

MVZinbaCopula::~MVZinbaCopula() {
  // free the memory of the correlation matrix
  free(this->rho);
}


void MVZinbaCopula::setSigma(double* sigma) {
  // compute the correlation matrix, so we do not have to do it everytime the density is evaluated 

  // the correlation coefficient in row I column J of the 
  // correlation matrixshould be stored in 
  // CORREL( J + ((I-2)*(I-1))/2 ), for J < I. 

  // in the R package the upper triangle is used, so we do this as well 

  // cout<<"alloc cor matrix: "<<this->p * (this->p - 1) / 2<<"\n";
  int k = 0;
  for (int j=1; j<this->p; j++) {
    for (int i=0; i<j; i++) {
      // cout<<"sigma["<<i<<", "<<j<<"] = "<<this->sigma[i + j * this->p]<<"\n";
      // cout<<"sigma["<<i<<", "<<i<<"] = "<<this->sigma[i + i * this->p]<<"\n";
      // cout<<"sigma["<<j<<", "<<j<<"] = "<<this->sigma[j + j * this->p]<<"\n";
      // cout<<"rho["<<i<<", "<<j<<"] = "<<this->sigma[i + j * this->p] / sqrt(this->sigma[i + i * this->p] * this->sigma[j + j * this->p])<<"\n";
      // cout<<"rho["<<k<<"] = ...\n";
      // sigma is stored by columns 
      rho[k++] = this->sigma[i + j * this->p] / sqrt(this->sigma[i + i * this->p] * this->sigma[j + j * this->p]);
      //k++;
    }
  }
}


double MVZinbaCopula::density(int t) {
  double* x = &(this->O[this->p * t]);
  return(density(x));
}

double MVZinbaCopula::density(double* x) {
 
#ifdef DEBUG
  int print_condition = (t >= 2352 && t <= 2357 || t < 5);
  if (print_condition) {
    cout<<"MVZinbaCopula::density("<<t<<") x = c(";
    for (int i=0; i<this->p-1; i++) {
      cout<<x[i]<<", ";
    }
    cout<<x[this->p - 1]<<") ";
  }
#endif

  // check if the density for this pair has been computed already 
  // precomputed values are stored in a hashmap 
//  map<int, map<int, double>*>::iterator it_x;
//  map<int, double>::iterator it_y;
//  it_x = precomputed->find((int)x);
//  // if not found the iterator is equal to end
//  if (it_x != precomputed->end()) { 
//    // the hashmap values are stored in the field second (keys in first)
//    it_y = it_x->second->find((int)y);
//    if (it_y != it_x->second->end()) {
//      return(it_y->second);
//    }
//  }
//
  // compute for a new pair 

  // upper limits for the integration 
  // first transform to p-values (uniform on [0, 1]) 
  // then transform to z-scores and store them in a double pointer 
  double* upper = (double*)calloc(this->p, sizeof(double));
  for (int i=0; i<this->p; i++) {
    double p = this->zinba[i]->CDF(&x[i]);
    // the last two arguments to qnorm are lower tail / logp
    upper[i] = qnorm(p, 0, 1, 1, 0) / sqrt(sigma[i + i * this->p]);
  }

  // we also need the lower limits for the integration 
  double* lower = (double*)calloc(this->p, sizeof(double));
  for (int i=0; i<this->p; i++) {
    double x_minus_one = x[i] - 1;
    double lp = this->zinba[i]->CDF(&x_minus_one);    
    lower[i] = qnorm(lp, 0, 1, 1, 0) / sqrt(sigma[i + i * this->p]);
  }

  /* evaluate the multivariate gaussian CDF using the following fortran routine
     (it is copied from the mvtnorm package: mvt.f)

      SUBROUTINE MVTDST( N, NU, LOWER, UPPER, INFIN, CORREL, DELTA, 
                         MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM )       

     A subroutine for computing non-central multivariate t probabilities.
     This subroutine uses an algorithm (QRSVN) described in the paper
     "Comparison of Methods for the Computation of Multivariate 
         t-Probabilities", by Alan Genz and Frank Bretz
         J. Comp. Graph. Stat. 11 (2002), pp. 950-971.

          Alan Genz 
          Department of Mathematics
          Washington State University 
          Pullman, WA 99164-3113
          Email : AlanGenz@wsu.edu

       Original source available from
       http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f

       This is version 7/10 with better support for 100 < dimension < 1000

     Parameters

     N      INTEGER, the number of variables.    
     NU     INTEGER, the number of degrees of freedom.
            If NU < 1, then an MVN probability is computed.
     LOWER  DOUBLE PRECISION, array of lower integration limits.
     UPPER  DOUBLE PRECISION, array of upper integration limits.
     INFIN  INTEGER, array of integration limits flags:
             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
     CORREL DOUBLE PRECISION, array of correlation coefficients; 
            the correlation coefficient in row I column J of the 
            correlation matrixshould be stored in 
               CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
            The correlation matrix must be positive semi-definite.
     DELTA  DOUBLE PRECISION, array of non-centrality parameters.
     MAXPTS INTEGER, maximum number of function values allowed. This 
            parameter can be used to limit the time. A sensible 
            strategy is to start with MAXPTS = 1000*N, and then
            increase MAXPTS if ERROR is too large.
     ABSEPS DOUBLE PRECISION absolute error tolerance.
     RELEPS DOUBLE PRECISION relative error tolerance.
     ERROR  DOUBLE PRECISION estimated absolute error, 
            with 99% confidence level.
     VALUE  DOUBLE PRECISION estimated value for the integral
     INFORM INTEGER, termination status parameter:
            if INFORM = 0, normal completion with ERROR < EPS;
            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
                           function vaules used; increase MAXPTS to 
                           decrease ERROR;
            if INFORM = 2, N > 1000 or N < 1.
            if INFORM = 3, correlation matrix not positive semi-definite.
  */

  int N = this->p;
  int df = 0;

  // set the infinity indicators 
  int* inf = (int*)calloc(this->p, sizeof(int));
  for (int i=0; i<this->p; i++) {
    if (x[i] == 0) {
      inf[i] = 0;
    } else {
      inf[i] = 2;
    }
  }

  
  double* delta = (double*)calloc(this->p, sizeof(double));

  // these are the defaults from the mvtnorm R package 
  int maxfuncalls = 25000;
  double releps = 0;
  double abseps = 0.001;
  double error = 0;
  int info = 0;
  double value = 0;

  // cout<<"call to mvtdst with upper = c("<<upper[0]<<", "<<upper[1]<<") ";
  // cout<<"lower = c("<<lower[0]<<", "<<lower[1]<<")\n";
  // cout<<"inf = c("<<inf[0]<<", "<<inf[1]<<")\n";
  

  F77_CALL(mvtdst)(&N, &df, lower, upper, inf, this->rho, delta, &maxfuncalls, &abseps, &releps, &error, &value, &info);
  
#ifdef DEBUG
  if (print_condition) {
    cout<<" value: "<<value<<" info: "<<info<<" error: "<<error;
  }
#endif

  if (ISNAN(value)) {
    cout<<"NAN value after mvtdst: "<<value<<" info: "<<info<<" error: "<<error<<"\n";
    cout<<"x = ";
    print_vector(x, this->p);
    cout<<"\ncall to mvtdst with\nupper = ";
    print_vector(upper, this->p);
    cout<<"\nlower = ";
    print_vector(lower, this->p);
    cout<<"\ninf = ";
    print_ivector(inf, this->p);
    cout<<"\n";
  }
  
  // there are problems with zero density values
  // if the integral gets very small, we use a rieman approximation
  // the area of the integration interval is one, so the using the normal
  // density at the upper limits is a lower bound (assuming we are in the
  // right tail of the distribution)
  if (value == 0 || ISNAN(value)) {
    value = mvn->density(upper);
  }

  // most likely one of the z values is Inf
  if (ISNAN(value)) {
    value = 0;
  }


  // we also observed some numerical instabilities where the value gets
  // negative even though integration limits should yield positive results
  if (value < 0) {
    // not sure what to do, just make it positive or approximate?
    value = fabs(value);
  }

#ifdef DEBUG
  if (print_condition) {
    cout<<" final value: "<<value<<"\n";
  }
#endif

  // free the auxillary variables 
  free(lower);
  free(upper);
  free(inf);
  free(delta);
  
  // store the results for later lookup 
//  if (it_x == precomputed->end()) { 
//    map<int, double>* ymap = new map<int, double>;
//    precomputed->insert(pair<int, map<int, double>*>((int)x, ymap));
//    it_x = precomputed->find((int)x);
//  }
//  it_x->second->insert(pair<int, double>((int)y, value));
//
  return(value);
}

double MVZinbaCopula::logdensity(int i) {
  return(0);
}

double MVZinbaCopula::logdensity(double* x) {
  return(0);
}

double MVZinbaCopula::CDF(double* x) {
  cout<<"MVZinbaCopula::CDF not implemented yet!\n";
  return(0);
}

double MVZinbaCopula::logCDF(double* x) {
  return(log(CDF(x)));
}


void MVZinbaCopula::initialize(double* weight, int T) {
  
}

void MVZinbaCopula::put(std::ostream& os) {
  os<<"MVZinbaCopula(size=c(";
  for (int i=0; i<p - 1; i++) {
    os<<size[i]<<", ";
  }
  os<<size[p - 1]<<"), mu=c(";
  for (int i=0; i<p - 1; i++) {
    os<<mu[i]<<", ";
  }
  os<<mu[p - 1]<<"), beta=c(";
  for (int i=0; i<p - 1; i++) {
    os<<beta[i]<<", ";
  }
  os<<beta[p - 1]<<"), sigma=c(";
  for (int i=0; i<p * p - 1; i++) {
    os<<sigma[i]<<", ";
  }
  os<<sigma[p * p - 1]<<"));\n";
}

void MVZinbaCopula::copy(Density* other) {
  
}

void MVZinbaCopula::update(double* weight, int T) {
  
};

DensityName MVZinbaCopula::getType() {
  return(Other);
}

double MVZinbaCopula::getMean() {
  return(0);
}

/* for testing of the mv copula density */
extern "C" {
  void R_dmvzinbacopula(double* O, int* N, int* p, double* size, double* mu, double* beta, double* sigma, double* density) {
    MVZinbaCopula* zc = new MVZinbaCopula(O, size, mu, beta, sigma, *p);

    cout<<*zc<<"\n";

    for (int i=0; i<*N; i++) {
      density[i] = zc->density(i);
    }
  }
}


MVCopulaApproximation::MVCopulaApproximation(double* observations, vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant) {
  this->O = observations;
  // these are the marginal distributions (we need their CDF function)
  this->marginals = marginals;
  this->Nmod = this->marginals.size();
  this->cor_matrix_inv = cor_matrix_inv;
  this->cor_matrix_determinant = cor_matrix_determinant;
}


double MVCopulaApproximation::logdensity(double* Ot) {
  double sum, exponent, exponentTemp;
  sum = 0.0;
  for(int imod=0;imod<Nmod;imod++){
    //this should return the value of the univariate density for the marginal, at position t, for that modification
    sum += this->marginals[imod]->logdensity(&Ot[imod]);
  }

  double* z = (double*)calloc(Nmod, sizeof(double));
  for(int imod=0; imod<Nmod; imod++) {
    // not sure that using long double here helps if CDF is only double..
    double uniform = this->marginals[imod]->CDF(&(Ot[imod]));//& because we are pointing to the element imod of Ot
    z[imod] = qnorm(uniform, 0, 1, 1, 0); // call to inverse erf in boost 
    //TODO: change last argument of qnorm to 1, and then provide loguniform, built the function logCDF for both zinb and NB
  }

/*cout<<endl;
cout<<"WHAT WOULD BE THE INVERSE OF THE CORR MATRIX"<<endl;
for(int imod=0;imod<Nmod;imod++){
for(int jmod=0;jmod<Nmod;jmod++){
cout<<this->cor_matrix_inv[imod * Nmod + jmod]<<"\t";
};cout<<endl;
};
cout<<endl;*/

  exponent = 0.0;
  for(int imod=0;imod<Nmod;imod++){
    exponentTemp = 0.0;
    for(int jmod=0;jmod<Nmod;jmod++){
      if(imod==jmod) exponentTemp+=z[jmod] * (this->cor_matrix_inv[imod * Nmod + jmod]- 1 );
      else exponentTemp += z[jmod] * this->cor_matrix_inv[imod * Nmod + jmod];
/*cout<<"(imod,jmod)=("<<imod<<","<<jmod<<"), exponentTemp = ";
if(imod==jmod) cout<<z[jmod]<<" * ( "<<this->cor_matrix_inv[imod * Nmod + jmod]- 1<<" )";
else cout<<z[jmod]<<" * "<<this->cor_matrix_inv[imod * Nmod + jmod];
cout<<" = "<<exponentTemp<<endl;*/
    }
    exponent+=exponentTemp * z[imod];
//cout<<"after loop jmod, exponent="<<exponent<<endl;
  }
  free(z);
//cout<<"exponent="<<exponent<<endl;
  return log(1 / sqrt(this->cor_matrix_determinant)) - 0.5 * exponent + sum;
}


double MVCopulaApproximation::density(double* Ot) {
  return(exp(this->logdensity(Ot)));
}

double MVCopulaApproximation::density(int t) {
  double* x = &(this->O[t * this->Nmod]);
  return(this->density(x));
}

double MVCopulaApproximation::logdensity(int t) {
  double* x = &(this->O[t * this->Nmod]);
  return(this->logdensity(x));
}

void MVCopulaApproximation::update(double* weight, int T) {
    //return 0;

}

double MVCopulaApproximation::getMean() {
  return(0);
}

DensityName MVCopulaApproximation::getType() {
  return(Other);
}
void MVCopulaApproximation::initialize(double *a, int b){
    //return 0;
}

void MVCopulaApproximation::put(ostream& s)
{
   // return 0;
}
void MVCopulaApproximation::copy(Density* d)
{
   // return 0;
}
double MVCopulaApproximation::CDF(double *d)
{
    return 0;
}
double MVCopulaApproximation::logCDF(double* x) {
  return(log(CDF(x)));
}


// Maria's univariate distributions
ZiNB::ZiNB() {
}

ZiNB::ZiNB(double* observations, double r, double p, double w)
{
    this->O = observations;//now should be only the reads, not the model
    this->p = p;
    this->r = r;
    this->w = w;
    this->cdf = new vector<double>();
    this->max_x = -1;//0;

    cout<<"Warning: CDF function not tested\n";
};

double ZiNB::density(double* Ot)
{
  double negbino;
  double lGammaR,lGammaRplusX,lxfactorial;
  lGammaRplusX=lgamma(this->r + *Ot);
  lGammaR=lgamma(this->r);
  lxfactorial=0;
  if (*Ot > 1){//if reads=0, log(xfactorial)=log(1)=0; if reads=1, log(xfactorial)=log(1)=0
		for(int j=2;j<=*Ot;j++){
			lxfactorial+=log(j);
		};
	};
	negbino=exp(lGammaRplusX-lGammaR-lxfactorial)*pow(this->p,this->r)*pow(1-this->p,*Ot);
	if(*Ot==0) return this->w+(1-this->w)*negbino;
	else return (1-this->w)*negbino;
};

double ZiNB::logdensity(double* Ot)
{
    double negbino;
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaRplusX=lgamma(this->r+*Ot);
	lGammaR=lgamma(this->r);
	lxfactorial=0;
	if(*Ot>1){//if reads=0, log(xfactorial)=log(1)=0; if reads=1, log(xfactorial)=log(1)=0
		for(int j=2;j<=*Ot;j++){
			lxfactorial+=log(j);
		};
	};
    if(*Ot==0)
    {
//if(isnan(log( this->w + (1-this->w) * exp(lGammaRplusX-lGammaR-lxfactorial) * pow(this->p,this->r)))){
//cout<<" reads="<<*Ot<<", logd=log("<<this->w<<" + "<<(1-this->w)<<" * exp("<<lGammaRplusX<<"-"<<lGammaR<<"-"<<lxfactorial<<") * "<<pow(this->p,this->r)<<"))"<<endl;
//};
        return( log( this->w + (1-this->w) * exp(lGammaRplusX-lGammaR-lxfactorial) * pow(this->p,this->r)));
    };
    if(*Ot>0)
    {
//if(isnan(log(1-this->w) + lGammaRplusX-lGammaR-lxfactorial + this->r * log(this->p) + *Ot * log(1-this->p) )){
//cout<<" reads="<<*Ot<<", logd= log("<<1-this->w<<") + "<<lGammaRplusX<<"-"<<lGammaR<<"-"<<lxfactorial<<" + "<<this->r <<"*"<< log(this->p)<<" + "<<*Ot <<"*log("<<1-this->p<<") )"<<endl;
//};
        return( log(1-this->w) + lGammaRplusX-lGammaR-lxfactorial + this->r * log(this->p) + *Ot * log(1-this->p) );
    };
//	negbino=exp(lGammaRplusX-lGammaR-lxfactorial)*pow(this->p,this->r)*pow(1-this->p,*Ot);
//	if(*Ot==0) return log(this->w+(1-this->w)*negbino);
//	else return log((1-this->w)*negbino);
    return(0);
};


double ZiNB::density(int t) {
  double* x = &(this->O[t]);
  return(this->density(x));
}

double ZiNB::logdensity(int t) {
  double* x = &(this->O[t]);
//if(t==2626) cout<<"t="<<t<<", reads="<<x<<endl;
  return(this->logdensity(x));
}

double ZiNB::CDF(double* Ot)
{
	int k = (int) *Ot;
	if (k < 0) return(0);
	// check if we have already computed the cdf this far
	if (k <= max_x) {
		double value = this->cdf->at(k);
		return(value);
	} else {
           // removed dependency on gsl and use pbeta instead

 	   //double lGamma1plusRplusX,lHyper,lGammaR,lGamma2plusX,lppowert,lppowerr;
           double current = 0;
           // lGammaR=lgamma(this->r);
           // lppowerr=this->r * log(this->p);

	   for (int i=max_x+1; i<=k; i++)
           {
    	   
	     
                // lGamma1plusRplusX=lgamma(1 + this->r + i);
                // lGamma2plusX=lgamma(2 + i);
		// 
                // lHyper=log(gsl_sf_hyperg_2F1(1, 1 + this->r + i, 2 + i, 1-this->p));
		// 
                // lppowert=(1 + i) * log(1-this->p);
	        // current = 1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX );
		// if(current == 1) current = this->cdf->at(i-1);

	        current = (1 - w) * pbeta(p, r, (double)i + 1, 1, 0);
		if (i == 0) {
		  current = w + current; // zero inflation
		}
            	cdf->push_back(current);
	   }
           max_x = k;
	   return(current);
	}
/*		double lGamma1plusRplusX,lHyper,lGammaR,lGamma2plusX,lppowert,lppowerr;
		long double x;
		lGamma1plusRplusX=lgamma(1+this->r+*Ot);
		lGammaR=lgamma(this->r);
		lGamma2plusX=lgamma(2+*Ot);
		lHyper=log(gsl_sf_hyperg_2F1(1, 1+this->r+*Ot, 2+*Ot, 1-this->p));
		lppowert=(1+*Ot)*log(1-this->p);
		lppowerr=this->r*log(this->p);
		x = 1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX );
		return ((double) x);
*/
}

double ZiNB::logCDF(double* x) {
  return(log(CDF(x)));
}


//TODO: in the loghmm, need to call posterior() function at each loop, before update, to get the weights
void ZiNB::update(double* weight, int T)
{
    int k, kmax, jN;
	double Fp, Fr, dr, where, phere, rhere, tempp, tempr, temp2p, temp2r;
	double dFpdp, dFrdr;
	double DigammaR, DigammaRplusX, DigammaRplusDRplusX, DigammaRplusDR;
	const double eps=1e-5;
	kmax=20;
	dr=0.00001;
	where=this->w;
	phere=this->p;
	rhere=this->r;
//	cout<<"weight="<<where<<", p="<<phere<<", r="<<rhere<<endl;
/*cout<<"weight"<<endl;
for(k=0;k<10;k++){
cout<<weight[k]<<endl;
};*/
	///////////////////////////////////////////////////////////////////////
    /////   NEWTON METHOD
    ///////////////////////////////////////////////////////////////////////
	for(k=1;k<kmax;k++){
        Fp=Fr=dFpdp=dFrdr=0.0;
	DigammaR = digamma(rhere); //boost::math::digamma<>(rhere);
	DigammaRplusDR = digamma(rhere + dr); //boost::math::digamma<>(rhere+dr);
		tempp=where+(1-where)*pow(phere,this->r);
		tempr=where+(1-where)*pow(this->p,rhere);
		temp2p=pow(where+(1-where)*pow(phere,this->r),2);
		temp2r=pow(where+(1-where)*pow(this->p,rhere),2);
		/*for(int ti=0;ti<T;ti++){
			sum=0;
			for(jN=0;jN<N;jN++){
				tempvec[jN]=this->model->logalpha[ti][jN]+this->model->logbeta[ti][jN];
			};
			temp=Max(tempvec, N);
			for(int jN=0;jN<N;jN++){
				sum+=exp(this->model->logalpha[ti][jN]+this->model->logbeta[ti][jN]-temp);
			};
			Tauti[ti]=exp(this->model->logalpha[ti][iN]+this->model->logbeta[ti][iN]-temp)/sum;//if proba=0 then Tauti=0
		};*/
		for(int ti=0;ti<T;ti++){
		  DigammaRplusX = digamma(rhere+this->O[ti]); //boost::math::digamma<>(rhere+this->O[ti]);
		  DigammaRplusDRplusX = digamma(rhere+dr+this->O[ti]); //boost::math::digamma<>(rhere+dr+this->O[ti]);
			if(this->O[ti]==0){
				Fp+=weight[ti]/tempp*(1-where)*pow(phere,this->r-1)*this->r;
                Fr+=weight[ti]/tempr*(1-where)*pow(this->p,rhere)*log(this->p);
                dFpdp+=-weight[ti]*(where-1)*pow(phere,this->r-2)*this->r*((where-1)*pow(phere,this->r)+where*(this->r-1))/temp2p;
                dFrdr+=-weight[ti]*(where-1)*where*pow(this->p,rhere)*pow(log(this->p),2)/temp2r;
			};
			if(this->O[ti]!=0){
				Fp+=weight[ti]*(this->r/phere-this->O[ti]/(1-phere));
                Fr+=weight[ti]*(log(this->p)-DigammaR+DigammaRplusX);
                dFpdp+=-weight[ti]*(this->r/pow(phere,2)+this->O[ti]/pow(phere-1,2));
                dFrdr+=weight[ti]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
			};
		};
		if(fabs(Fp)+fabs(Fr)<eps) break;
		if(Fp/dFpdp<phere) phere=phere-Fp/dFpdp;
		if(Fp/dFpdp>phere) phere=phere/2.0;
		if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
		if(Fr/dFrdr>rhere) rhere=rhere/2.0;
		if(k%50==0){
			cout<<"NewtonRaphson: k= "<<k<<", (Fp,Fr)= "<<Fp<<", "<<Fr;
			cout<<", (dp,dr)="<<-Fp/dFpdp<<", "<<-Fr/dFrdr;
			cout<<", (p,r)="<<phere<<", "<<rhere<<endl;
		};
	};
	this->p=phere;
	this->r=rhere;
//cout<<"------------------ UPDATES: r= "<<this->r<<" , p= "<<this->p<<" , w= "<<this->w<<endl;
}

void ZiNB::initialize(double* weight, int T) {
  cout<<"ZiNB::update not implemented yet!\n";
}

double ZiNB::getMean()
{
	return (1-this->w)*this->r*(1-this->p)/this->p;
}

DensityName ZiNB::getType() {
  return(Z);
}

void ZiNB::put(ostream &os) {
  os << "ZiNB(p = " << this->p << ", r = " << this->r << ", w = " << this->w <<")";
}

void ZiNB::setR (double newR) {
  this->r = newR;
}

double ZiNB::getR() {
  return(this->r);
}

void ZiNB::setP (double newP) {
  this->p = newP;
}

double ZiNB::getP() {
  return(this->p);
}

void ZiNB::setW (double newW) {
    this->w = newW;
}

double ZiNB::getW() {
  return(this->w);
}

void ZiNB::copy(Density* other) {
  cout<<"ZiNB::copy";
  ZiNB* o = (ZiNB*)other;
  this->p = o->p;
  this->r = o->r;
  this->O = o->O;
  this->w = o->w;
  cout<<" OK!\n";
  
}

NegativeBinomial::NegativeBinomial() {}

NegativeBinomial::NegativeBinomial(double* observations, double r, double p){//Hmm *model, double r, double p) {
    this->O = observations;
    this->r = r;
    this->p = p;
    this->cdf = new vector<double>();
    this->max_x = -1;//0;
    
    cout<<"Warning: CDF function not tested\n";
};

double NegativeBinomial::density(double* Ot) {
  double negbino;
  double lGammaR,lGammaRplusX,lxfactorial;
  /*TODO: maybe a good idea to pre-compute the logarithm of the factorial for the largest number of reads*/
  lGammaRplusX=lgamma(this->r+*Ot);
  lGammaR=lgamma(this->r);
  lxfactorial=0;
  if(*Ot>1){//if reads=0, log(xfactorial)=log(1)=0; if reads=1, log(xfactorial)=log(1)=0
    for(int j=2;j<=*Ot;j++){
      lxfactorial+=log(j);
    };
  };
  negbino=exp(lGammaRplusX-lGammaR-lxfactorial)*pow(this->p,this->r)*pow(1-this->p,*Ot);
  return negbino;
};

double NegativeBinomial::logdensity(double* Ot){//Ot) {
  double negbino;
  double lGammaR,lGammaRplusX,lxfactorial;
  lGammaRplusX=lgamma(this->r+*Ot);
  lGammaR=lgamma(this->r);
  lxfactorial=0;
  if(*Ot>1){//if reads=0, log(xfactorial)=log(1)=0; if reads=1, log(xfactorial)=log(1)=0
    for(int j=2;j<=*Ot;j++){
      lxfactorial+=log(j);
    }
  }
  negbino=lGammaRplusX-lGammaR-lxfactorial + this->r * log(this->p) + *Ot * log(1-this->p);
//if(isnan(negbino)) 
//cout<<"negative binomial, reads="<<Ot<<", density:"<<endl;
//cout<<"logd="<<lGammaRplusX<<"-"<<lGammaR<<"-"<<lxfactorial<<" +"<< this->r <<"*"<< log(this->p) <<"+"<< *Ot <<"*"<< log(1-this->p)<<endl;
//cout<<"logd="<<negbino<<endl;
//};
  return negbino;
};

double NegativeBinomial::density(int t) {
  double* x = &(this->O[t]);
  return(this->density(x));
}

double NegativeBinomial::logdensity(int t) {
  double* x = &(this->O[t]);
//if(t==2626) cout<<endl<<"reads="<<x<<endl;
  return(this->logdensity(x));
//cout<<"what I return is      this->logdensity(x) = "<<this->logdensity(x)<<endl;
}

double NegativeBinomial::CDF(double* Ot)
{
  int k = (int) (*Ot);
  if (k < 0) return(0);
  // check if we have already computed the cdf this far
  if (k <= max_x) {////////////////PROBLEM HERE IF k=0 at the first bin and max_x = 0!!!!!!!!!!
    double value = this->cdf->at(k);//cdf[k]
    return(value);
  } else {
        // removed dependency on gsl and use pbeta instead

	// double lGamma1plusRplusX,lHyper,lGammaR,lGamma2plusX,lppowert,lppowerr;
    	// lGammaR=lgamma(this->r);
        // lppowerr=this->r * log(this->p);
	double current = 0;
    	for (int i=max_x+1; i<=k; i++) /////////////PROBLEM HERE!!!!!!!!!! I START AT 1, NOT AT ZERO
	{
	    // lGamma1plusRplusX=lgamma(1 + this->r + i);
	    // lGamma2plusX=lgamma(2 + i);
	    // lHyper=log(gsl_sf_hyperg_2F1(1, 1 + this->r + i, 2 + i, 1-this->p));
	    // lppowert=(1+i) * log(1-this->p);
	    // current = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX );
	    // if(current == 1) {cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);};

	    current = pbeta(p, r, (double)i + 1, 1, 0);
      	    cdf->push_back((double) current);
	}
    	max_x = k;
    	return(current);
  }
/*    double lGamma1plusRplusX,lHyper,lGammaR,lGamma2plusX,lppowert,lppowerr;
    long double x;
    lGamma1plusRplusX=lgamma(1+this->r+k);
    lGammaR=lgamma(this->r);
    lGamma2plusX=lgamma(2+k);
    lHyper=log(gsl_sf_hyperg_2F1(1, 1+this->r+k, 2+k, 1-this->p));
    lppowert=(1+k)*log(1-this->p);
    lppowerr=this->r*log(this->p);
    x = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX );
    return ((double) x);
*/
}

double NegativeBinomial::logCDF(double* x) {
  return(log(CDF(x)));
}


//TODO: in the loghmm, need to call posterior() function at each loop, before update, to get the weights
void NegativeBinomial::update(double* weight, int T) {
    int ti, jN, k, kmax;
    double eps = 1e-4;
    double numerator, denominator, rhere,dr,Fr,dFrdr,DigammaR,DigammaRplusDR,DigammaRplusX,DigammaRplusDRplusX;
    ///////////////// UPDATE p
    numerator=denominator=0.0;
    for(ti=0;ti<T;ti++){
        numerator+=weight[ti]*this->r;
        denominator+=weight[ti]*(this->r+this->O[ti]);
    };
    this->p = numerator/denominator;//////CAUTION! THE UPDATE OF r IS NOW DONE WITH THIS UPDATED p
    ///////////////// UPDATE r
    ///////////////   NEWTON METHOD
    ///////////////////////////////////////////////////////////////////////
	rhere=this->r;
	dr=0.00001;
	kmax=20;
    for(k=1;k<kmax;k++){
        Fr=dFrdr=0.0;
        DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
        DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
        for(int ti=0;ti<T;ti++){
	  DigammaRplusX = digamma(rhere+this->O[ti]); //boost::math::digamma<>(rhere+this->O[ti]);
	  DigammaRplusDRplusX = digamma(rhere+dr+this->O[ti]); // boost::math::digamma<>(rhere+dr+this->O[ti]);
            if(this->O[ti]==0){
                Fr+=weight[ti]*log(this->p);
                //dFrdr+=0;
            };
            if(this->O[ti]!=0){
                Fr+=weight[ti]*(log(this->p)-DigammaR+DigammaRplusX);
                dFrdr+=weight[ti]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
            };
        };
        if(fabs(Fr)<eps) break;
        if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
        if(Fr/dFrdr>rhere) rhere=rhere/2.0;
        if(k%50==0){
            cout<<"NewtonRaphson: k= "<<k<<", Fr= "<<Fr<<", dr="<<-Fr/dFrdr<<", r="<<rhere<<endl;
        };
    };
    this->r=rhere;
//cout<<"------------------ UPDATES: r= "<<this->r<<" , p= "<<this->p<<endl;
	//p[iN]= numerator/denominator;
};

void NegativeBinomial::initialize(double* weight, int T) {
  cout<<"ZiNB::update not implemented yet!\n";
}

double NegativeBinomial::getMean()
{
  return this->r*(1-this->p)/this->p;
}

DensityName NegativeBinomial::getType() {
  return(NB);
}

void NegativeBinomial::put(ostream &os) {
  os << "NegativeBinomial(p = " << this->p << ", r = " << this->r << ")";
}



void NegativeBinomial::setR (double newR) {
  this->r = newR;
}

double NegativeBinomial::getR() {
  return(this->r);
}

void NegativeBinomial::setP (double newP) {
  this->p = newP;
}

double NegativeBinomial::getP() {
  return(this->p);
}

void NegativeBinomial::copy(Density* other) {
  cout<<"NegativeBinomial::copy";
  NegativeBinomial* o = (NegativeBinomial*)other;
  this->p = o->p;
  this->r = o->r;
  this->O = o->O;
  cout<<" OK!\n";
}
