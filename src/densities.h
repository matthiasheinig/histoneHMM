#ifndef DENSITIES_H
#define DENSITIES_H


#include <Rcpp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
#include <fstream>
using std::ifstream;
using std::ostream;
#include <cstdlib> 

// for some strage reason, this header has to be loaded after Rcpp
// #include <Rinternals.h>

// #include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>

// #include <map>
#include "zeroin.h"


// Rcpp actually renames all rfunction to Rf_.. when using the namespace R
// we do not need to change all function names! Just define the std stuff above

#include <R.h>

using namespace R;


// for the multivariate normal CDF we need to interface with some fortran code
extern "C" { 
  void F77_NAME(mvtdst)(int* N, int* df, double* lower, double* upper, int* inf, double* cor, double* mu, int* maxfuncalls, double* abseps, double* releps, double* value, double* error, int* info);
}

enum DensityName {Other, Z, NB};

class Density {
 public:
  virtual double density(int t) = 0; // t is the index of the data point
  /* this is the density for the actual data point x 
     (it has to be a pointer to be general enough for the multivariate case) */
  virtual double density(double* x) = 0;
  // log densities
  virtual double logdensity(int t) = 0;
  virtual double logdensity(double* x) = 0;
  // this is the cumulative density
  virtual double CDF(double* x) = 0;
  // log of the CDF
  virtual double logCDF(double* x) = 0;
  // this is for the EM / Baum-Welch algorithm
  virtual void update(double* weight, int T) = 0; 
  // this is to compute initial parameter estimates
  virtual void initialize(double* weight, int T) = 0;
  // this is the function to print to the output stream
  virtual void put(ostream& os) = 0;
  // this overrides the << operator for the printing
friend ostream& operator <<(ostream &os, Density &obj);
  virtual void copy(Density* other) = 0;
  // get the type of distribution
  virtual DensityName getType() = 0;
  // get the mean
  virtual double getMean() = 0;
  /* observations (if it is multivariate use column major encoding 
     O[ndim x nobs] */
  double* O;
  int n; // number of observations
};  

class Zero : public Density {
  //  friend ostream& operator <<(ostream &os, const LogNormal &obj);
 public:
  Zero();
  Zero(double* observations);
  double density(int t);
  double density(double* x);
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
 private:
  double* O; // observations
};

class LogNormal : public Density {
  //  friend ostream& operator <<(ostream &os, const LogNormal &obj);
 public:
  LogNormal();
  LogNormal(double* observations, int n, double mu, double sigma);
  double density(int t);
  double density(double* x);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  double mu;
  double sigma;
  DensityName getType();
  double getMean();
 private:
  double* O; // observations
};


/* we also want bivariate distributions for the differential states */
class BivariateLogNormal : public Density {
  //  friend ostream& operator <<(ostream &os, const BivariateLogNormal &obj);
 public:
  BivariateLogNormal();
  BivariateLogNormal(double* observations, int n, double* mu, double* sigma, double rho);
  double density(int t);
  double density(double* x);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
 private:
  double* O; // observations
  double* mu;
  double* sigma;
  double rho;
};


class BivariateNormal : public Density {
 public:
  BivariateNormal();
  BivariateNormal(double* observations, int n, double* mu, double* sigma, double rho);
  BivariateNormal(double* observations, double mu_x, double sigma_x, double mu_y, double sigma_y, double rho);
  double density(int t);
  double density(double* x);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
 private:
  double* O; // observations
  double mu_x;
  double sigma_x;
  double mu_y;
  double sigma_y;
  double rho;
};


class MultivariateNormal : public Density {
 public:
  MultivariateNormal(double* observations, double* mu, double* sigma, int p);
  ~MultivariateNormal();
  double density(int t);
  double density(double* x);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  void setSigma(double* sigma);
  DensityName getType();
  double getMean();
 private:
  double* O; // observations
  double* mu;
  double* sigma;
  double* sigma_inv;
  double sigma_det;
  int p;
};

void matrixMultiply(int transposeA, int transposeB, 
		    int nrowA, int ncolA, int nrowB, int ncolB, 
		    double *a, double *b, double *c);

/* zero inflated negative binomial distribution */

class Zinba : public Density {
  //  friend ostream& operator <<(ostream &os, const Density &obj);
 public:
  Zinba();
  Zinba(double* observations, int n, double size, double mu, double beta);
  Zinba(double size, double mu, double beta);
  ~Zinba();
  double density(int t);
  double density(double* x);
//  double logLikelihood();
//  double Zinba_objective(int n, double *par, void *ex);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  static double density(double x, double size, double mu, double beta);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void initialize(double* weight);
  void put(std::ostream& os);
  void copy(Density* other);
  double logLikelihood();
  void update(double* weight, int T);
  double size;
  double mu;
  double beta;
  DensityName getType();
  double getMean();
 private:
  double* O;
  int n; // size of O
  vector<double>* logcdf;
  int max_x;
  double* weight; // this stores the weights for updates (only temporarily)
};


// objective function for the numerical optimization
double Zinba_objective(int n, double *par, void *ex);

/* zero inflated negative binomial distribution */
class Negbinom : public Density {
  //  friend ostream& operator <<(ostream &os, const Density &obj);
 public:
  Negbinom();
  Negbinom(double* observations, int n, double size, double mu);
  Negbinom(double size, double mu);
  ~Negbinom();
  double density(int t);
  double density(double* x);
  static double density(double x, double size, double mu);
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize();
  void initialize(double* weight, int T);
  DensityName getType();
  double getMean();
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  void maria_update(double* weight, int T);
  double logLikelihood();
  double size;
  double mu;

 private:
  double* O;
  int n; // size of O
  vector<double>* logcdf;
  int max_x;
  double* weight; // this stores the weights for updates (only temporarily)
  static double logLikelihoodPartialDerivative(double size, void* extra);
};

/* bivariate zero inflated negative binomial

   try if we can slice the data along one axis and fit the conditional distr.
   which is also a negative binomial
  
   we use a formulation of the form P(X,Y) = P(X|Y) P(Y)
   where P(Y) is a zero inflated negative binomial and
   P(X|Y) is also a zero inflated negative binomial.
   for the estimation we create equavilly sized bins
   since we know that in a bivariate NB the conditional distribution 
   P(X|Y) ~ NB(size + y, p) we use a linear model to compute the size

*/


class BivariateZinba : public Density {
 public:
  BivariateZinba();
  BivariateZinba(double* observations, int n, double marginal_size, double marginal_mu, double marginal_beta, double* size_coef, double* mu_coef, double beta0);
  double density(int t);
  double density(double* x);
  // log densities
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
  double marginal_size;
  double marginal_mu;
  double marginal_beta;
  double beta0;
  double* size_coef;
  double* mu_coef;
};


/* here is a more elegant way using a gaussian copula on top of the zero 
   inflated negative binomials for the marginal distributions */


class ZinbaCopula : public Density {
 public:
  ZinbaCopula(double* O, double size_x, double mu_x, double beta_x, double size_y, double mu_y, double beta_y, double sigma_x, double sigma_y, double rho);
  ZinbaCopula(double* O, double size_x, double mu_x, double beta_x, double size_y, double mu_y, double beta_y, double sigma_x, double sigma_y, double rho, double copula_mu_x, double copula_mu_y);
  double density(int t);
  double density(double* x);
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
  double size_x;
  double mu_x;
  double beta_x;
  double size_y;
  double mu_y;
  double beta_y;
  double sigma_x;
  double sigma_y;
  double rho;
  double copula_mu_x;
  double copula_mu_y;

 private:
  // these are the marginal distributions (we need their CDF function)
  Zinba* zinba_x;
  Zinba* zinba_y;
  // we also keep a bivariate normal density to get lower bounds when the 
  // numerical accuracy of the normal CDF is getting too low
  BivariateNormal* bvn;
  // map<int, map<int, double>*>* precomputed;

};




/* here is a more elegant way using a gaussian copula on top of the zero 
   inflated negative binomials for the marginal distributions */
class MVZinbaCopula : public Density {
 public:
  MVZinbaCopula(double* O, double* size, double* mu, double* beta, double* sigma, int p);
  ~MVZinbaCopula();
  double density(int t);
  double density(double* x);
  double logdensity(int t);
  double logdensity(double* x);
  double CDF(double* x);
  double logCDF(double* x);
  void initialize(double* weight, int T);
  void put(std::ostream& os);
  void copy(Density* other);
  void update(double* weight, int T);
  DensityName getType();
  double getMean();
  void setSigma(double* sigma);
  double* size;
  double* mu;
  double* beta;
  double* sigma;

 private:
  // these are the marginal distributions (we need their CDF function)
  vector<Zinba*> zinba;
  // we also keep a multivariate normal density to get lower bounds when the 
  // numerical accuracy of the normal CDF is getting too low
  MultivariateNormal* mvn;
  double* mvn_mu;
  int p; // dimensions
  double* rho; // correlation matrix

};


class MVCopulaApproximation : public Density {
 public:
  MVCopulaApproximation(double* observations, vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant);
  ~MVCopulaApproximation();
  double density(double* Ot);
  double logdensity(double* Ot);
  double density(int t);
  double logdensity(int t);  
  void update(double* weight, int T);
  double getMean();
  DensityName getType();
    double CDF(double* t);
    double logCDF(double* x);
    void copy(Density* other);
    void initialize(double* weight, int T);
    void put(std::ostream &os);
 private:
  double* O; // observations
  vector<Density*> marginals;
  double* cor_matrix_inv;
  double cor_matrix_determinant;
  int Nmod;
};


// Maria's univariate densities
class ZiNB : public Density {
public:
  ZiNB();
  ZiNB(double* observations, double r, double p, double w);
  ~ZiNB();//TODO if something is allocated
  double density(double* Ot);
  double logdensity(double* Ot);
  double density(int t);
  double logdensity(int t);  
  double CDF(double* t);
  double logCDF(double* x);
  void update(double* weight, int T);
  void initialize(double* weight, int T);
  double getMean();
  DensityName getType();
  void put(std::ostream &os);
  void setR (double newR);
  double getR();
  void setP (double newP);
  double getP();
    void setW (double newW);
  double getW();
  void copy(Density* other);
 private:
  double* O; // observations
  double p;
  double r;
  double w;
  int max_x;
  vector<double>* cdf;
};


class NegativeBinomial : public Density {
public:
  NegativeBinomial();
  NegativeBinomial(double* observations, double r, double p);
  ~NegativeBinomial();//TODO if something is allocated
  double density(double* Ot);
  double logdensity(double* Ot);
  double density(int t);
  double logdensity(int t);  
  double CDF(double* x);
  double logCDF(double* x);
  void update(double* weight, int T);
  void initialize(double* weight, int T);
  double getMean();
  DensityName getType();
  void put(std::ostream &os);
  void setR (double newR);
  double getR();
  void setP (double newP);
  double getP();
  void copy(Density* other);
 private:
  double* O; // observations
  double p;
  double r;
  int max_x;
  vector<double>* cdf;
};




#endif
