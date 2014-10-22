#ifndef _MEANSHIFT_
#define _MEANSHIFT_

#include <math.h>

typedef struct {
  double * points;
  double * thresh;
  int iterations;
  double * start;
  double * final;

} MeanShiftRep;

double * profileMd(double *a, double *b, int n, double *x, double h);
double * profile1d(double *xi, int n, double x, double h);
double * meanshift(double *a, double *b, int n, double *x, double h);
double euclidian_norm(double *x, int d);
MeanShiftRep meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter);

#endif
