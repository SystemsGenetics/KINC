#ifndef _SWILK_
#define _SWILK_

#include <math.h>
#include <gsl/gsl_statistics.h>
#include <stdio.h>
#include <string.h>
#include "stats.h"
#include "../vector.h"

double poly(const double *cc, int nord, double x);
void swilk(double *vector, int n, double *w, double *pw, int *ifault);

#endif
