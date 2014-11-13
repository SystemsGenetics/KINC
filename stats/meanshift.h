#ifndef _MEANSHIFT_
#define _MEANSHIFT_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stats.h"

typedef struct {
  // Maps the points in the 2D vector with the cluster number they belong to.
  // The first element corresponds to the first point, the second element to
  // the second point, etc.
  int * cluster_label;
  // Indicates the total number of clusters
  int num_clusters;
  // An array that indicates the size of each cluster. The first element
  // is the size of the cluster labeled 1, the second element is the size of
  // the cluster labeled 2, etc.
  int * sizes;
  // Organizes the points by cluster.  The first element corresponds to
  // the cluster labeled 1, the second to the cluster labeled 2, etc.  Each
  // element is itself an array of points from the 2D vector.
  double *** clusters;

} MeanShiftClusters;

typedef struct {
  double ** points;
  double * thresh;
  int iterations;
  double * start;
  double * final;

} MeanShiftRep;

MeanShiftClusters meanshift2D(double* x, double * y, int n, double h);
MeanShiftRep meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter);

// Helper functions
double * profileMd(double *a, double *b, int n, double *x, double h);
double * profile1d(double *xi, int n, double x, double h);
double * meanshift_base(double *a, double *b, int n, double *x, double h);
double euclidian_norm(double *x, int d);
int minimal_dist(double **x, int n, double *y);
double * distance_vector(double **x, int n, double *y);

#endif
