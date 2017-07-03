
#ifndef HMM_H
#define HMM_H

#include "utility.h"
#include "densities.h"



/* we store all info related to the model in this struct */
class HMM  {

public:
  HMM();
  HMM(double* O, int T, int N);
  HMM(double* O, int T, int N, int Nmod);
  ~HMM();

  double* O; /* sequence of observations */
  int T; /* length of observed sequence */
  int N; /* number of states */
  int Nmod; /* number of dimensions (modifications/marks) */
  double** A; /* : matrix [N x N] of transition probabilities */
  double* proba; /* initial probabilities (length N) */
  vector<Density*> densityFunctions; /* density functions for each state */


  void initialize();
  double logLikelihood();
  void computeDensities();
  void forward();
  void backward();
  void estimateTransitions(int iterationMAX, double eps);
  void baumWelch(int iterationMAX, double eps);
  void viterbi(int* path, int recompute);
  void posterior(double** post, int recompute);

};


/* implementation using the scaling strategy */ 
class ScalingHMM  {

public:
  ScalingHMM();
  ScalingHMM(double* O, int T, int N);
  ScalingHMM(double* O, int T, int N, int Nmod);
  ~ScalingHMM();

  void initialize();
  double logLikelihood();

  double* O; /* sequence of observations */
  int T; /* length of observed sequence */
  int N; /* number of states */
  int Nmod; /* number of dimensions (histone modifications) */
  double** A; /* : matrix [N x N] of transition probabilities */
  double* proba; /* initial probabilities (length N) */
  vector<Density*> densityFunctions; /* density functions for each state */
  void computeDensities();
  void forward();
  void backward();
  void estimateTransitions(int iterationMAX, double eps);
  void baumWelch(int iterationMAX, double eps);
  void viterbi(int* path, int recompute);
  void posterior(double** post, int recompute);

private:
  double sumGamma(int begin,  int end, int i);
  double sumEta(int begin, int end, int i, int j);

  double** alpha; /* matrix [T x N] of forward probabilities */
  double** beta; /*  matrix [T x N] of backward probabilities */
  double** alphanonorm; /* matrix [T x N] of forward probabilities not scaled */
  double** betanonorm; /* matrix [T x N] of backward probabilities not scaled */
  double** densities; /* matrix [T x N] of precomputed density values */
};


/* implementation using the scaling strategy */
class LogHMM  {
    
public:
    LogHMM();
    LogHMM(double* O, int T, int N);
    LogHMM(double* O, int T, int N, int Nmod);
    
    ~LogHMM();
    
    void initialize();
    double logLikelihood();
    void computeDensities();
    double* O; /* sequence of observations */
    int T; /* length of observed sequence */
    int N; /* number of states */
    int Nmod; //number of modifications/marks
    double** A; /* : matrix [N x N] of transition probabilities */
    double* logproba; /* initial probabilities (length N) */
    vector<Density*> densityFunctions; /* density functions for each state */
    
    void forward();
    void backward();
    int estimateTransitions(int iterationMAX, double eps); //return the last iteration so that we know how many iterations it actually runs through
    int baumWelch(int iterationMAX, double eps); //return the last iteration so that we know how many iterations it actually runs through
    void viterbi(int* path, int recompute);
    void posterior(double** post, int recompute);
    
    bool zn2zz; //for the univariate HMM if it needs to change the density functions from (Z, NB) to (Z, Z),  in the construction zn2zz is set to FALSE by default, the value will be changed to TRUE in the baumWelch if necessary
    
private:
    double logSumGamma(int i);
    double logSumEta(int i, int j);
    
    double** logalpha; /* matrix [T x N] of forward probabilities */
    double** logbeta; /*  matrix [T x N] of backward probabilities */
    double** logdensities; /* matrix [T x N] of precomputed logdensity values */
};

#endif
