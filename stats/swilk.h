#ifndef _SWILK_
#define _SWILK_

#include <gsl/gsl_statistics_float.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "stats.h"
#include "../general/vector.h"
#include "../general/error.h"

float poly(const float *cc, int nord, float x);
void swilk(float *vector, int n, float *w, float *pw, int *ifault);

#endif
