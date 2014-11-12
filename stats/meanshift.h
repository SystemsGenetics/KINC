#ifndef _MEANSHIFT_
#define _MEANSHIFT_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stats.h"

typedef struct {
  double ** points;
  double * thresh;
  int iterations;
  double * start;
  double * final;

} MeanShift;

int * meanshift2D(double* x, double * y, int n, double h);

double * profileMd(double *a, double *b, int n, double *x, double h);
double * profile1d(double *xi, int n, double x, double h);
double * meanshift_base(double *a, double *b, int n, double *x, double h);

double euclidian_norm(double *x, int d);
MeanShift meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter);

int minimal_dist(double **x, int n, double *y);
double * distance_vector(double **x, int n, double *y);

#endif
