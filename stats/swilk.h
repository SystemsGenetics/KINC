#ifndef _SWILK_
#define _SWILK_

#include <gsl/gsl_statistics.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "stats.h"
#include "../general/vector.h"
#include "../general/error.h"

double poly(const double *cc, int nord, double x);
void swilk(double *vector, int n, double *w, double *pw, int *ifault);

#endif
