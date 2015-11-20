#ifndef _SFRANCIA_
#define _SFRANCIA_

#include <math.h>
#include <gsl/gsl_statistics.h>
#include <stdio.h>
#include <string.h>
#include "stats.h"
#include "../general/vector.h"
#include "../general/error.h"

void sfrancia(double *vector, int n, double *w, double *pw, int *ifault);

#endif
