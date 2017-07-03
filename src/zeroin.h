#ifndef ZEROIN_H
#define ZEROIN_H

#include "math.h"
#include "float.h"

#include <iostream>

/* An estimate to the root	*/
double zeroin(double ax, double bx, double (*f)(double x, void* extra), double tol, void* extra);
/* double ax;				 Left border | of the range	*/
/* double bx;  			         Right border| the root is seeked*/
/* double (*f)(double x, void* extra)    Function under investigation	*/
/* double tol;				 Acceptable tolerance		*/

#endif
