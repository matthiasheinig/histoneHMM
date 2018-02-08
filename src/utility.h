#ifndef UTILITY_H
#define UTILITY_H


#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
using std::ifstream;
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <float.h>
#include <time.h>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
//MA
#include <exception>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <map>
#include <sstream>

//Matthias's density.h
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
// #include <Rinternals.h>
#include <R_ext/Utils.h> //from Matthias's hmm.h
//to use GSL library
// #include <gsl/gsl_sf_hyperg.h>//needed for the hypergeometric function

//#include <boost/locale.hpp>
//#include </home/p265768/boost_1_51_0/boost/math/distributions/negative_binomial.hpp>
//#include </home/p265768/boost_1_51_0/boost/math/special_functions/digamma.hpp>

//for MinhAnh
// #include <boost/algorithm/string.hpp>
// #include <boost/math/distributions/negative_binomial.hpp>
// #include <boost/math/special_functions/digamma.hpp>
// #include <boost/math/distributions/gamma.hpp>#
// #include <boost/math/special_functions/beta.hpp>

//for Maria
//#include </home/p261756/EURATRANS/boost_1_51_0/boost/math/special_functions/beta.hpp>
//#include </home/p261756/EURATRANS/boost_1_51_0/boost/algorithm/string.hpp>
//#include </home/p261756/EURATRANS/boost_1_51_0/boost/math/special_functions/digamma.hpp>
//#include </home/p261756/EURATRANS/boost_1_51_0/boost/math/distributions/gamma.hpp>

// 

// /**============ CONSTANT VARIABLES ============================================*/
const double pi = M_PI;
//enum DensityName {NB, Z}; //NB = Negative Binomial, Z = ZiNBa
// 
// /**========== GENERAL UTILITIES ===================================================*/
// /**
// 	print error message then exit program
// */
// void outError(const char *error);
// /**
//   	print error message then exit program
// */
// void outError(string error);
// /**
// 	print double error messages then exit the program
// */
// void outError(const char *error, const char *msg);
// 
// /**
// 	print double error messages then exit program
// */
// void outError(const char *error, string msg);
// /**
// 	convert string to integer, with error checking
// 	@param str original string
// 	@return the integer value
// */
// int convert_int(const char *str) throw (string);
// /**
// 	convert string to double, with error checking
// 	@param str original string
// 	@return the double
// */
// double convert_double(const char *str) throw (string);
// /**
//  * convert int to string
//  * @param int
//  * @return string
//  */
// string convertIntToString(int number);
// 
// /**
//  * convert all character in a string to lower chacacters
//  * @param str the input string
//  * @RETURN the string with all characters are "lowered"
//  */
// string convertToLower(const string str);
// /**
//  * search a string in a vector of strings, return the stringID if found or -1 if not found
//  * @param str the input string
//  * @param str_vector the given string vector
//  * @RETURN stringID or -1
//  */
// //int searchString(const string str, const StringVec& str_vector);
// 
// /**
//  * return the ordinary round of a non-integer number
//  */
// int ordinaryRounding(double a);
// 
// /**
//  * FOR THE HMM FROM MATTHIAS
//  **/

/* helpers for memory management */
double** allocDoubleMatrix(int rows, int cols);
void freeDoubleMatrix(double** matrix, int rows);

// 
// int** allocIntMatrix(int rows, int cols);
// void freeIntMatrix(int** matrix, int rows);
// 
 double*** alloc3Ddouble(int dim1, int dim2, int dim3);
 void free3Ddouble(double*** array, int dim1, int dim2);
// 
 bool** allocBoolMatrix(int rows, int cols);
 void freeBoolMatrix(bool** matrix, int rows);
// 
// //double sumP(int begin, int end, double** x, double* norm, int N, int T);
// /**
//  * Helpers to find maximum, an index of the maximum value
//  */
//for a double array
double Max(double *a, int N); //return the maximum value
int argMax(double *a, const int N); //return an index of the maximum (the first one if tight happens)
//for an integer array
int intMax(int *a, int N); //return the maximum value
int argIntMax(int *a, const int N); ////return an index of the maximum (the first one if tight happens)
//for a double matrix
double MaxMatrix(double**, int N, int M);//return the maximum value
//for an integer matrix
int MaxIntMatrix(int**, int N, int M);//return the maximum value
double MaxDoubleMatrix(double**, int N, int M);

// void readUniDensityParam(const char *infile, double **densparam);
// void readUniCall(const char *infile, bool **enrich, const int markID);
// 


#endif // UTILITY_H
