#include "meanshift.h"

/**
 * Calculates the Mean Shift Clusters (MCL) of two vectors.
 *
 * This function was adapted from the R function ms() which is
 * available in the royston package at:
 * http://cran.r-project.org/web/packages/LPCM/index.html
 *
 * Unlike the R code, this function has been adapted to only deal with
 * two vectors rather than a matrix, and has fewer arguments.
 *
 * @param double * a
 * @param double * b
 * @param int n
 *   The size of a and b.
 * @param h
 *   The bandwidth: the percentage of the range to use for clustering.
 */
void meanshift2D(double* a, double * b, int n, double h, double thr, int iter) {
  int d = 2; // the number of columns (i.e. only two vectors).

  // First calculate the size of the range for each vector
  int i = 0;
  int si[2];
  int a_min = INFINITY;
  int b_min = INFINITY;
  int a_max = -INFINITY;
  int b_max = -INFINITY;

  for (i = 0; i < n; i++) {

    if (a_min < a[n]) {
      a_min = a[n];
    }
    if (a_max > a[n]) {
      a_max = a[n];
    }

    if (b_min < b[n]) {
      b_min = b[n];
    }
    if (b_max > b[n]) {
      b_max = b[n];
    }
  }
  si[0] = a_max - a_min;
  si[1] = b_max - b_min;

  double finals[n][d] = 0;
  int ncluster = 0;
  MeanShiftRep temp_msr;

  // Iterate through the vectors
  for (i = 0; i < n; i++) {
    temp_msr = meanshift_rep(a, b, n, h, 1e-8, iter);
  }

}

/**
 * Mean shift iterative function (until convergence ...)
 */
MeanShiftRep meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter) {
  int d = 2;
  int s = 0;
  int j = 0;
  double *x0 = x;
  double **M;
  double th[iter] = 0;
  double *m;
  double mx[d];
  MeanShiftRep msr;

  M = (double **) malloc(sizeof(double *) * iter);
  for (j = 0; j < iter; j++) {
    m = meanshift(a, b, n, x, h);
    M[j] = m;
    mx[0] = m[0] - x[0];
    mx[1] = m[1] - x[1];
    th[j] = euclidian_norm(&mx) / euclidian_norm(x);
    if (th[j] < thresh){
      s = j;
      break;
    }
    x = m;
  }


  msr.points = M;
  msr.iterations = s;
  msr.start = x0;
  msr.final = m;
  msr.thresh = th;

  return msr;
}

/**
 * Euclidian norm
 */
double euclidian_norm(double *x, int d) {
  int i;
  double sum = 0;
  for (i = 0; i < d; i++) {
    sum += pow(x[i], 2);
  }
  return sum;
}

/**
 * Mean shift base function
 *
 * @param double *a
 * @param double *b
 * @param double n
 *   The size of a and b.
 * @param double *x
 *   An array of only two values containing representing a "row" from
 *   the 2D matrix of a and b.
 * @param double h
 *   The bandwidth: the percentage of the range to use for clustering.
 *
 * @return double *
 *   An double array with two elements.
 */
double * meanshift(double *a, double *b, int n, double *x, double h) {
  int d = 2;
  int i;
  double * g;
  double * ms =  (double *) malloc(sizeof(double) * 2);
  double sum_g = 0;
  double sum_ag = 0;
  double sum_bg = 0;
  memset(ms, 0, n);

  g = gd(a, b, n, x, h);
  for (i = 0; i < n; i++) {
    sum_g += g[i];
    sum_ag += a[i] * g[i];
    sum_bg += b[i] * g[i];
  }
  ms[0] = sum_ag / sum_g;
  ms[1] = sum_bg / sum_g;
  free(g);
  return ms;
}

/**
 *  Calculate the 1-d profile.
 *
 *  @param double * xi
 *  @param int n
 *  @param double x
 *  @param double h
 *
 *  @return double *
 *    An array of length n containing the 1 dimensional profile
 */
double * profile1d(double *xi, int n, double x, double h) {
  double *k1 = (double *) malloc(sizeof(double) * n);
  int j;
  for (j = 0; j < n; j++) {
    k1[j] = 1/2 * exp(-1/2 * ((x-xi[j]) / h)^2);
  }
  return k1;
}

/**
 * Calculate the Multi-d profile.
 *
 * @param double *a
 * @param double *b
 * @param double n
 *   The size of a and b.
 * @param double *x
 *   An array of only two values containing representing a "row" from
 *   the 2D matrix of a and b.
 * @param double h
 *   The bandwidth: the percentage of the range to use for clustering.
 *
 * @return double *
 *    An array of length n containing the multi-dimensional profile
 */
double * profileMd(double *a, double *b, int n, double *x, double h) {
  int d = 2; // only 2 dimensions in our matrix (a & b)
  int i, j;
  double *pa, *pb;
  double *k = (double *) malloc(sizeof(double) * n);

  // Intialize the k array to contain all 1's
  memset(k, '1', n);

  // Get the 1D profiles for both a and b. We use the corresponding value
  // of x to pass in to the profile1d function (a first, b second).
  pa = profile1d(a, n, x[0], h);
  pb = profile1d(b, n, x[1], h);

  // Calculate the multidimensional profile
  for (i = 0; i < n; i++) {
    k[i] = k[i] * pa[i];
    k[i] = k[i] * pb[i];
  }

  free(pa);
  free(pb);

  return k;
}
