#ifndef _MEANSHIFT_
#define _MEANSHIFT_

/**
 * The code for Mean Shift Clustering and bandwidth selection was adapted
 * from the LPCM R page found here:
 * http://artax.karlin.mff.cuni.cz/r-help/library/LPCM/html/ms.html
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../stats/stats.h"
#include "../vector.h"
#include "PairWiseClusters.h"

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
  float *** clusters;
  // The center points of each cluster.
  float ** centers;

  // The scaled input data vectors
  float * a;
  float * b;

  float scaled_by[2];

} MeanShiftClusters;

typedef struct {
  float ** points;
  float * thresh;
  int iterations;
  float * start;
  float * final;
} MeanShiftRep;

MeanShiftClusters * meanshift2D(float* x, float * y, int n, float h);
MeanShiftRep * meanshift_rep(float* a, float * b, int n, float * x, float h,
    float thresh, int iter);

void free_msr(MeanShiftRep * msr);
void free_msc(MeanShiftClusters * msc);

// Helper functions
float * profileMd(float *a, float *b, int n, float *x, float h);
float * profile1d(float *xi, int n, float x, float h);
float * meanshift_base(float *a, float *b, int n, float *x, float h);
float euclidian_norm(float *x, int d);
float minimal_dist(float **x, int n, float *y);
float * distance_vector(float **x, int n, float *y);
float * meanshift_coverage2D(float *s, float *t, int n);
float coverage_raw(float * a, float *b, int n, float ** centers, int nc, float tau);

#endif
