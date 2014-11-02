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
void meanshift2D(double* a, double * b, int n, double h) {
  int iter = 200;
  double thr = 0.0001;

  // First calculate the size of the range for each vector
  int i, j;
  double si[2];
  double a_min = INFINITY;
  double b_min = INFINITY;
  double a_max = -INFINITY;
  double b_max = -INFINITY;

  for (i = 0; i < n; i++) {

    if (a_min > a[i]) {
      a_min = a[i];
    }
    if (a_max < a[i]) {
      a_max = a[i];
    }

    if (b_min > b[i]) {
      b_min = b[i];
    }
    if (b_max < b[i]) {
      b_max = b[i];
    }
  }
  si[0] = a_max - a_min;
  si[1] = b_max - b_min;

  double *finals[n];
  int ncluster = 0;
  double *savecluster[n];
  int cluster_label[n];
  int closest_label[n];
  double cluster_dist[n];
  MeanShift temp_ms;
  double min_dist = 0;
  double which_min = 0;
  double x[2];


  // Iterate through the elements of a & b (i.e the rows of the 2D matrix)
  for (i = 0; i < n; i++) {

    // Create the point, x.
    x[0] = a[i];
    x[1] = b[i];

    temp_ms = meanshift_rep(a, b, n, x, h, 1e-8, iter);
    finals[i] = temp_ms.final;
    cluster_dist[n] = 0;

    // If we have a cluster then calculate the distance of this point to the
    // cluter and find the one  with the minimum distance.  This will be the
    // cluster to which the point, x, belongs
    if (ncluster >= 1) {
      min_dist = INFINITY;
      for (j = 0; j < ncluster; j++) {
        double t[2];
        t[0] = savecluster[j][0] - finals[i][0];
        t[1] = savecluster[j][1] - finals[i][1];
        cluster_dist[j] = euclidian_norm(t, 2) / euclidian_norm(savecluster[j], 2);
        if (min_dist < cluster_dist[j]) {
          min_dist = cluster_dist[j];
          which_min = j;
        }
      }
    }
    // If we have no clusters or the minimum distance is greater than the
    // threshold then perform the following.
    if (ncluster == 0 || min_dist > thr) {
      ncluster = ncluster + 1;
      savecluster[ncluster] = finals[i];
      cluster_label[i] = ncluster;
    }
    // If we have a cluster and the minimum distance is less than the threshold
    // then add a label to the point
    else {
      cluster_label[i] = which_min;
    }
  }

  for (i = 0; i < n; i++){
    // Create the point, x.
    x[0] = a[i];
    x[1] = b[i];
    int md = minimal_dist(savecluster, ncluster, x);
    closest_label[i] = md;
  }
}

/**
 * Mean shift iterative function (until convergence ...)
 */
MeanShift meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter) {
  int d = 2;
  int s = 0;
  int j = 0;
  double *x0 = x;
  double **M;
  double th[iter];
  double *m;
  double mx[d];
  MeanShift msr;

  M = (double **) malloc(sizeof(double *) * iter);
  for (j = 0; j < iter; j++) {
    m = meanshift(a, b, n, x, h);
    M[j] = m;
    mx[0] = m[0] - x[0];
    mx[1] = m[1] - x[1];
    th[j] = euclidian_norm(mx, 2) / euclidian_norm(x, 2);
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
 * Computes the minimal distance between a vector and a set of vectors.
 *
 * @param double **x
 *   A 2D array of doubles, or a list of points.
 * @param int n
 *   The length of x
 * @param double *y
 *
 * @return double *
 *   A two dimensional array where the first element is the minimum
 *   distance and the second is the position in the array where the
 *   minimum distance is found.
 */
int minimal_dist(double **x, int n, double *y) {
  int d = n;
  int i;
  double * dv = distance_vector(x, n, y);
  double s[n];
  double min_s = INFINITY;
  int min_i = 0;

  for (i = 0; i < n; i++) {
    s[n] = sqrt(d) * dv[i];
    if (s[n] < min_s) {
      min_s = s[n];
      min_i = i;
    }
  }
  free(dv);

  // Return the minimum distance and the position.
//  double *out = (double *) malloc(sizeof(double) * 2);
//  out[0] = min_s;
//  out[1] = min_i;
  return min_i;
}

/**
 * Computes all distances between a set of vectors and another vector.
 * (up to the constant sqrt(d))
 */
double * distance_vector(double **x, int n, double *y) {
  int i;
  int d = 2;

  // Reference the two vectors in x as a and b.
  double *a = x[0];
  double *b = x[n];

  // The original R code looks like the following.
  // out <- sqrt(as.vector(rowMeans(X^2) + mean(y^2) - 2 * X %*% y/n)))
  // the %*% means matrix multiplication

  // Calculate the squared mean of each row in our ab Matrix, and the squared
  // mean of y.
  double x2means[n];
  double y2mean = 0;
  for (i = 0; i < n; i++) {
    x2means[n] = (pow(a[i], 2) + pow(b[i], 2)) / 2;
  }
  y2mean = (pow(y[0], 2) + pow(y[1],2)) / d;

  // Get a vector by  performing matrix multiplication between ab and y.
  // Multiply the value by 2 and divide by n.
  double Xy[n];
  for (i = 0; i < n; i++) {
    Xy[i] = 2 * (a[i] * y[0] + b[i] * y[1]) / n;
  }

  // Finish up the distance vector.
  double * out = (double *) malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    out[i] = sqrt(x2means[i] + y2mean - Xy[i]);
    if (isnan(out[i])) {
      out[i] = 0;
    }
  }

  return out;
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
  int i;
  double * g;
  double * ms =  (double *) malloc(sizeof(double) * 2);
  double sum_g = 0;
  double sum_ag = 0;
  double sum_bg = 0;
  memset(ms, 0, n);

  g = profileMd(a, b, n, x, h);
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
    k1[j] = 1/2 * exp(-1/2 * pow((x-xi[j]) / h, 2));
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
  int i;
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
