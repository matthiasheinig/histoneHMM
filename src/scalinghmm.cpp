#include "densities.h"
#include "hmm.h"
#include "utility.h"

// #define DEBUG

// using namespace std;

// namespace histoneHMM {

// const double pi = 3.14159265358979323846264338327950288419716939937510;

/* helpers for memory management */
// double** allocDoubleMatrix(int rows, int cols) {
//   double** matrix = (double**) calloc(rows, sizeof(double*));
//   int i;
//   for (i=0; i<rows; i++) {
//     matrix[i] = (double*) calloc(cols, sizeof(double));
//   }
//   return(matrix);
// }
// 
// void freeDoubleMatrix(double** matrix, int rows) {
//   int i;
//   for (i=0; i<rows; i++) {
//     free(matrix[i]);
//   }
//   free(matrix);
// }





/* initialize the model */
ScalingHMM::ScalingHMM(double* O, int T, int N) {
  this->O = O;
  this->T = T;
  this->N = N;
  this->A = allocDoubleMatrix(N, N);
  this->alpha = allocDoubleMatrix(T, N);
  this->beta = allocDoubleMatrix(T, N);
  this->alphanonorm = allocDoubleMatrix(T, N);
  this->betanonorm = allocDoubleMatrix(T, N);
  this->densities = allocDoubleMatrix(T, N);
  this->proba = (double*) calloc(N, sizeof(double));

  /* we set initial probabilities to uniform and make self transitions more
     likely than transitions to other states */
  double self = 0.9;
  double other = (1.0 - self) / (N - 1.0);
  for (int i=0; i<N; i++) {
    this->proba[i] = 1.0 / N;
    cout<<"proba["<<i<<"] = "<<this->proba[i]<<"\n";
    for (int j=0; j<N; j++) {
      if (i == j) {
	this->A[i][j] = self;
      } else {
	this->A[i][j] = other;
      }
      cout<<"A["<<i<<"]["<<j<<"] = "<<this->A[i][j]<<"\n";
    }
  }
}

ScalingHMM::ScalingHMM(double* O, int T, int N, int Nmod) {
  ScalingHMM(O, T, N);
  this->Nmod = Nmod;
}

/* clean up a model (do not touch the input data) */
ScalingHMM::~ScalingHMM() {
  freeDoubleMatrix(this->A, this->N);
  freeDoubleMatrix(this->alpha, this->T);
  freeDoubleMatrix(this->beta, this->T);
  freeDoubleMatrix(this->alphanonorm, this->T);
  freeDoubleMatrix(this->betanonorm, this->T);
  freeDoubleMatrix(this->densities, this->T);
  free(this->proba);
}

// just a dummy function
void ScalingHMM::initialize() {};

void ScalingHMM::computeDensities() {};

/* log likelihood */

/* x matrix [T x N] */
double sumP(int begin, int end, double** x, int N, int T){
  double s, sum;
  int q, i;
  s = 0.0;
  for (q=begin; q<end; q++) {
    sum = 0.0;
    for (i=0; i<N; i++) {
      sum += x[q][i];
    }
    if (sum == 0) {
      /* if likelihood is numerically zero we set it to the smallest float */
      sum = DBL_MIN;
      cout<<"\nsumP zero: "<<q<<" set to "<<sum<<" log(1/sum) = "<<log(1.0 / sum)<<"\n";
    }
    s += log(1.0 / sum);
  }
  return s;
}

/* helper functions */
double ScalingHMM::sumGamma(int begin,  int end, int i){
  double s, sum;
  int q, k;
  s = 0.0;
  for (q=begin;q<end;q++){
    sum = 0.0;
    for (k=0; k<this->N; k++) {
      sum += this->alpha[q][k] * this->beta[q][k];
    }
    s += this->alpha[q][i] * this->beta[q][i] / sum;
  }
  return s;
}

double ScalingHMM::sumEta(int begin, int end, int i, int j) {
  double s;
  int q, r, z;
  double SUM;
  s = 0.0;
  for(q=begin; q<end; q++) {
    SUM = 0.0;
    for (r=0; r<this->N; r++) {
      for (z=0; z<this->N; z++) {
	SUM += this->alpha[q][r] * this->A[r][z] * this->densities[q + 1][z] * this->beta[q + 1][z];
      }
    }
    s += this->alpha[q][i] * this->A[i][j] * this->densities[q + 1][j] * this->beta[q + 1][j] / SUM;
  }
  return s;
}


void ScalingHMM::forward() {

  /* printf("forward\n"); */

  double normalization, sum;
  int i, j, t;

  /* CALCULATE AND STORE THE ALPHAS */
  /* first row */
  for(i=0; i<this->N; i++) {
    this->alphanonorm[0][i] = this->proba[i] * this->densities[0][i];
  }
  normalization = 0.0;
  for(i=0; i<this->N; i++) {
    normalization += this->alphanonorm[0][i];
  }
  for(i=0; i<this->N; i++) {
    /* noramlized */
    this->alpha[0][i] = this->alphanonorm[0][i] / normalization;
  }
  /* recursion */
  for(t=1; t<this->T; t++) {
    for(i=0; i<this->N; i++) {
      sum=0.0;
      for(j=0; j<this->N; j++) {
	sum += this->alpha[t-1][j] * this->A[j][i];
      }
      this->alphanonorm[t][i] = sum * this->densities[t][i];
    }
    normalization = 0.0;
    for(i=0; i<this->N; i++) {
      normalization += this->alphanonorm[t][i];
    }
    for(i=0; i<this->N; i++) {
      /* noramlized */
      this->alpha[t][i] = this->alphanonorm[t][i] / normalization;
      /* if (t < 10) printf("%f ", this->alpha[t][i]); */
    }
    /* if (t < 10) printf("\n"); */
  }
}

void ScalingHMM::backward() {

  /* printf("backward\n"); */

  double normalization, sum;
  int i, j, t;
 
  /* CALCULATE AND STORE THE BETAS */
  /* last row */
  for(i=0; i<this->N; i++) {
    this->betanonorm[this->T - 1][i] = 1.0;
  }
  normalization = 0.0;
  for(i=0; i<this->N; i++) {
    normalization += this->alphanonorm[this->T - 1][i];
  }
  for(i=0; i<this->N; i++) {
    /* noramlized */
    this->beta[this->T - 1][i] = this->betanonorm[this->T - 1][i] / normalization;
  }
  /* recursion */
  for(t=this->T-2; t>=0; t--) {
    for(i=0; i<this->N; i++) {
      sum=0.0;
      for(j=0; j<this->N; j++) {
	sum += this->A[i][j] * this->densities[t + 1][j] * this->beta[t + 1][j];
      }
      this->betanonorm[t][i] = sum;
    }
    normalization = 0.0;
    for(i=0; i<this->N; i++) {
      normalization += this->alphanonorm[t][i];
    } 
    for(i=0; i<this->N; i++) {
      /* noramlized */
      this->beta[t][i] = this->betanonorm[t][i] / normalization;
      /* if (t < 10) printf("%f ", this->alpha[t][i]); */
    }
    /*  if (t < 10) printf("\n"); */
  }
}


void ScalingHMM::baumWelch(int iterationMAX, double eps) {
  
  /* arbitrary initialization was done in the constructor */

  double logPold = -INFINITY;
  double logPnew, smallnumber, sumDeno;
  double normalization, sum;
  /* useful number */
  double* SUMofGamma = (double*) calloc(this->N, sizeof(double));
  double** gamma = allocDoubleMatrix(this->N, this->T);
  double** Anew = allocDoubleMatrix(this->N, this->N);

  int i, j, k, t, iteration;

  for(iteration=1; iteration <= iterationMAX; iteration++) {
    R_CheckUserInterrupt(); // check interupt
    cout<<"iteration = "<<iteration;

    /* compute the density values */
    for(t=0; t<this->T; t++){
#ifdef DEBUG
      int print_condition = (t >= 2352 && t <= 2357 && iteration == 1);
      // int print_condition = (iteration == 1 && t < 10) {
      if (print_condition) {
	printf("density[%i,] = c(", t + 1);
      }
#endif
      for(i=0; i<this->N; i++){
	this->densities[t][i] = this->densityFunctions[i]->density(t);
#ifdef DEBUG
	if (print_condition) {
	  printf("%e", this->densities[t][i]);
	  if (i < this->N - 1) {
	    printf(", ");
	  }
	}
#endif
      }
#ifdef DEBUG
      if (print_condition) {
	printf(")\n");
      }
#endif
    }
    
    /* compute forward probabilities */
    forward();
    
    /* CALCULATE THE LOG LIKELIHOOD */
    logPnew = -sumP(0, this->T, this->alphanonorm, this->N, this->T);
    smallnumber = logPnew - logPold; /* fabs(logPnew - logPold); */
    
    cout<<" delta logP = " <<smallnumber<<"\n";

    /* check convergence */
    if(fabs(smallnumber) < eps){
      break;
    } else {
      logPold = logPnew;
    }
    
    /* compute backward probabilities */
    backward();

    /* updates */

    /* useful number */
    for(i=0; i<this->N; i++){
      SUMofGamma[i] = sumGamma(0, this->T, i);
    }
    
    /* update A */
    for(i=0; i<this->N; i++) {
      for(j=0; j<this->N; j++) {
	sumDeno = 0.0;
	for(k=0; k<this->N; k++) {
	  sumDeno += this->alpha[this->T - 1][k] * this->beta[this->T - 1][k];
	}
	Anew[i][j] = sumEta(0, this->T - 1, i, j) / (SUMofGamma[i] - this->alpha[this->T - 1][i] * this->beta[this->T - 1][i] / sumDeno);
      }
    }
    /* copy the new A to the model */
    for(i=0; i<this->N; i++) {
      for(j=0; j<this->N; j++) {
	this->A[i][j] = Anew[i][j];
      }
    }

    /* update initial probabilities */
    normalization = 0.0;
    for(i=0; i<this->N; i++) {
      normalization += this->alpha[0][i] * this->beta[0][i];
    }
    for(i=0; i<this->N; i++) {
      this->proba[i] = this->alpha[0][i] * this->beta[0][i] / normalization;
    }

    /* update the distributions in each state */
    cout<<"\nupdating state emissions:\n";

    /* compute the gamma (prob of being in state i at position t)
       and use this as weights for the update of the state's emissions */
    posterior(gamma, 0);
    for (i=0; i<this->N; i++) {
      cout<<i<<"\n";
      this->densityFunctions[i]->update(gamma[i], this->T);
    }
  } /* main loop end */

  /* free memory */
  free(SUMofGamma);
  freeDoubleMatrix(gamma, this->N);
  freeDoubleMatrix(Anew, this->N);
}



void ScalingHMM::estimateTransitions(int iterationMAX, double eps) {
  
  /* arbitrary initialization was done in the constructor */

  double logPold = -INFINITY;
  double logPnew, smallnumber, sumDeno;
  double normalization, sum;
  /* useful number */
  double* SUMofGamma = (double*) calloc(this->N, sizeof(double));
  double** Anew = allocDoubleMatrix(this->N, this->N);

  int i, j, k, t, iteration;

  /* compute the density values (in this function we can do this before the
     iterations since emission probabilities are fixed */
  for(t=0; t<this->T; t++){
#ifdef DEBUG
    int print_condition = (t >= 2352 && t <= 2357 && iteration == 1);
    // int print_condition = (iteration == 1 && t < 10) {
    if (print_condition) {
      printf("density[%i,] = c(", t + 1);
    }
#endif
    for(i=0; i<this->N; i++){
      this->densities[t][i] = this->densityFunctions[i]->density(t);
#ifdef DEBUG
      if (print_condition) {
	printf("%e", this->densities[t][i]);
	if (i < this->N - 1) {
	  printf(", ");
	}
      }
#endif
    }
#ifdef DEBUG
    if (print_condition) {
      printf(")\n");
    }
#endif
  }
 

  for(iteration=1; iteration <= iterationMAX; iteration++) {
    R_CheckUserInterrupt(); // check interupt
    cout<<"iteration = "<<iteration;
    
     
    /* compute forward probabilities */
    forward();
    
    /* CALCULATE THE LOG LIKELIHOOD */
    logPnew = -sumP(0, this->T, this->alphanonorm, this->N, this->T);
    smallnumber = logPnew - logPold; /* fabs(logPnew - logPold); */
    
    cout<<" delta logP = " <<smallnumber<<"\n";

    /* check convergence */
    if(fabs(smallnumber) < eps){
      break;
    } else {
      logPold = logPnew;
    }
    
    /* compute backward probabilities */
    backward();

    /* updates */

    /* useful number */
    for(i=0; i<this->N; i++){
      SUMofGamma[i] = sumGamma(0, this->T, i);
    }
    
    /* update A */
    for(i=0; i<this->N; i++) {
      for(j=0; j<this->N; j++) {
	sumDeno = 0.0;
	for(k=0; k<this->N; k++) {
	  sumDeno += this->alpha[this->T - 1][k] * this->beta[this->T - 1][k];
	}
	Anew[i][j] = sumEta(0, this->T - 1, i, j) / (SUMofGamma[i] - this->alpha[this->T - 1][i] * this->beta[this->T - 1][i] / sumDeno);
      }
    }
    /* copy the new A to the model */
    for(i=0; i<this->N; i++) {
      for(j=0; j<this->N; j++) {
	this->A[i][j] = Anew[i][j];
      }
    }

    /* update initial probabilities */
    normalization = 0.0;
    for(i=0; i<this->N; i++) {
      normalization += this->alpha[0][i] * this->beta[0][i];
    }
    for(i=0; i<this->N; i++) {
      this->proba[i] = this->alpha[0][i] * this->beta[0][i] / normalization;
    }
  } /* main loop end */

  /* free memory */
  free(SUMofGamma);
  freeDoubleMatrix(Anew, this->N);
}

/* for an initialized model compute the viterbi path */
void ScalingHMM::viterbi(int* path, int recompute) {
  if (recompute) {
    /* compute forward and backward variables */
    forward();
    backward();
  }
  
  /* do the actual viterbi */

}

/* for an initialized model compute the posterior probs 
   posterior is a matrix [T x N] */
void ScalingHMM::posterior(double** post, int recompute) {
  if (recompute) {
    /* compute forward and backward variables */
    forward();
    backward();
  }

  int i, t;
  double sum;

  for(t=0; t<this->T; t++) {
    sum = 0.0;
    for(i=0; i<this->N; i++) {
      sum += this->alpha[t][i] * this->beta[t][i];
    }
    for(i=0; i<this->N; i++) {
      post[i][t] = this->alpha[t][i] * this->beta[t][i] / sum;
      // if (t < 10) printf("(P=%f, a=%f, b=%f) ", post[t][i], this->alpha[t][i], this->beta[t][i]);
    }
    // if (t < 10) printf("\n");
  }  
}

/* compute the total loglikelihood of the data */
double ScalingHMM::logLikelihood() {
  return(-sumP(0, this->T, this->alphanonorm, this->N, this->T));
}


/* the functions that are called from R are defined here */
extern "C" {

  void R_hmm_posterior(double* O, int* T, int* N, double* mu, double* sigma, int* iterationMAX, double* eps, double* post, int* bivariate, double* A, double* proba, double* rho) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f, bivariate %i\n", *T, *N, *iterationMAX, *eps, *bivariate);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      // Density* d;
      if (*bivariate == 0) {
	printf("using lognormal\n");
	LogNormal* d = new LogNormal(O, *T, mu[i], sigma[i]);
	model->densityFunctions.push_back(d);
      } else {
	/* if the observations are bivariate the params come as two rows
	   of a column encoded matrix (for now we do not use correlation) */
	printf("using bivariate lognormal\n");
	BivariateLogNormal* d = new BivariateLogNormal(O, *T, &mu[i * 2], &sigma[i * 2], rho[i]);
	model->densityFunctions.push_back(d);
      }
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->estimateTransitions(*iterationMAX, *eps);
    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }

  void R_hmm_baumwelch(double* O, int* T, int* N, double* mu, double* sigma, int* iterationMAX, double* eps, double* post, int* bivariate, double* A, double* proba, double* rho) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f, bivariate %i\n", *T, *N, *iterationMAX, *eps, *bivariate);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      // Density* d;
      if (*bivariate == 0) {
	printf("using lognormal\n");
	LogNormal* d = new LogNormal(O, *T, mu[i], sigma[i]);
	model->densityFunctions.push_back(d);
      } else {
	/* if the observations are bivariate the params come as two rows
	   of a column encoded matrix (for now we do not use correlation) */
	printf("using bivariate lognormal\n");
	BivariateLogNormal* d = new BivariateLogNormal(O, *T, &mu[i * 2], &sigma[i * 2], rho[i]);
	model->densityFunctions.push_back(d);
      }
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->baumWelch(*iterationMAX, *eps);

    /* get back the parameters to R */
    for (i=0; i<model->N; i++) {
      // Density* d;
      if (*bivariate == 0) {
	mu[i] = ((LogNormal*)(model->densityFunctions[i]))->mu;
	sigma[i] = ((LogNormal*)(model->densityFunctions[i]))->sigma;
      } 
      // in the bivariate case we do not have to do anything because we passed
      // the pointers directly, so values should be up to date already
    }

    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }


void R_hmm_baumwelch_negbinom(double* O, int* T, int* N, double* mu, double* size, int* iterationMAX, double* eps, double* post, double* A, double* proba)
{
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      // check if a zero inflation component is included (encoded by mu = 0)
      if (mu[i] == 0) {
	Zero * d = new Zero(O);
	model->densityFunctions.push_back(d);
      } else {
	Negbinom * d = new Negbinom(O, *T, size[i], mu[i]);
	model->densityFunctions.push_back(d);
      }
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->baumWelch(*iterationMAX, *eps);

    /* get back the parameters to R */
    for (i=0; i<model->N; i++) {
      // Density* d;
      // only consider the non-zero distributions for copying
      if (mu[i] != 0) {
	mu[i] = ((Negbinom*)(model->densityFunctions[i]))->mu;
	size[i] = ((Negbinom*)(model->densityFunctions[i]))->size;
      }
    }

    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }

  /* run with zero inflated negative binomial distributions */
  void R_hmm_posterior_zinba(double* O, int* T, int* N, double* mu, double* size, double* beta, int* iterationMAX, double* eps, double* post, double* A, double* proba) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      printf("using zinba\n");
      Zinba* d = new Zinba(O, *T, size[i], mu[i], beta[i]);
      model->densityFunctions.push_back(d);
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->estimateTransitions(*iterationMAX, *eps);
    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }


  /* run with bivariate zero inflated negative binomial distributions */
  void R_hmm_posterior_bivariatezinba(double* O, int* T, int* N, double* marginal_mu, double* marginal_size, double* marginal_beta, double* size_coef, double* mu_coef, double* beta0, int* iterationMAX, double* eps, double* post, double* A, double* proba) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      printf("using bivariate zinba\n");
      BivariateZinba* d = new BivariateZinba(O, *T, marginal_size[i], marginal_mu[i], marginal_beta[i], &size_coef[i * 2], &mu_coef[i * 2], beta0[i]);
      model->densityFunctions.push_back(d);
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->estimateTransitions(*iterationMAX, *eps);
    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }


  /* run with copula version of zero inflated negative binomial distributions */
  void R_hmm_posterior_zinbacopula(double* O, int* T, int* N, double* size_x, double* mu_x, double* beta_x, double* size_y, double* mu_y, double* beta_y, double* sigma_x, double* sigma_y, double* rho, int* iterationMAX, double* eps, double* post, double* A, double* proba) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      printf("using zinba copula\n");
      ZinbaCopula* d = new ZinbaCopula(O, size_x[i], mu_x[i], beta_x[i], size_y[i], mu_y[i], beta_y[i], sigma_x[i], sigma_y[i], rho[i]);
      cout<<*d<<"\n";
      model->densityFunctions.push_back(d);
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->estimateTransitions(*iterationMAX, *eps);
    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }

  /* run with copula version of zero inflated negative binomial distributions */
  void R_hmm_posterior_mvzinbacopula(double* O, int* T, int* N, int* p, double* size, double* mu, double* beta, double* sigma, int* iterationMAX, double* eps, double* post, double* A, double* proba) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    /* printf("init..\n"); */
    ScalingHMM* model = new ScalingHMM(O, *T, *N);

    /* copy the distribution params */
    int i, j, t;
    for (i=0; i<model->N; i++) {
      printf("using multivariate zinba copula\n");
      MVZinbaCopula* d = new MVZinbaCopula(O, &size[i * *p], &mu[i * *p], &beta[i * *p], &sigma[i * *p * *p], *p);
      cout<<*d<<"\n";
      model->densityFunctions.push_back(d);
    }
    
    /* printf("OK\n"); */

    /* estimate the parameters */
    /* printf("estimate transitions..\n"); */
    model->estimateTransitions(*iterationMAX, *eps);
    /* printf("OK\n"); */

    /* compute the posteriors and save results directly to the R pointer */
    /* we do not need to recompute the forward and backward variables since 
       they have been computed in the estimation already */
    /* printf("compute posterior.."); */
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, 0);
    /* recode into column representation */
    for (t=0; t<model->T; t++) {
      for (i=0; i<model->N; i++) {
	post[t + i * model->T] = post_matrix[i][t];
      }
    }
    freeDoubleMatrix(post_matrix, model->N);
    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++) {
      proba[i] = model->proba[i];
      for (j=0; j<model->N; j++) {
	A[i + j * model->N] = model->A[i][j];
      }
    }
    /* printf("OK\n"); */
  }
}
