#ifndef _MEANSHIFT_
#define _MEANSHIFT_

#include <math.h>
#include "stats.h"

typedef struct {
  double * points;
  double * thresh;
  int iterations;
  double * start;
  double * final;

} MeanShift;

double * profileMd(double *a, double *b, int n, double *x, double h);
double * profile1d(double *xi, int n, double x, double h);
double * meanshift(double *a, double *b, int n, double *x, double h);
double euclidian_norm(double *x, int d);
MeanShift meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter);

double * min_dist(double *a, double *b, n, double *y);
double * distance_vector(double *a, double *b, n, double *y);

#endif
