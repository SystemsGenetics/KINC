#ifndef _SFRANCIA_
#define _SFRANCIA_

#include <math.h>
#include <gsl/gsl_statistics_float.h>
#include <stdio.h>
#include <string.h>
#include "stats.h"
#include "../general/vector.h"
#include "../general/error.h"

void sfrancia(float *vector, int n, float *w, float *pw, int *ifault);

#endif
