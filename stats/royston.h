#ifndef _ROYSTON_
#define _ROYSTON_

/*
 * The code form the rosyton.c was adapted from the R code in the
 * royston package: http://cran.r-project.org/web/packages/royston/
 */

#include <math.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_cdf.h>
#include "stats.h"
#include "kurtosis.h"
#include "sfrancia.h"
#include "swilk.h"
#include "../general/error.h"

float roystonH(float* a, float * b, int n, float *pcc);

#endif
