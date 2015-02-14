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
 * @param double * s
 * @param double * t
 * @param int n
 *   The size of a and b.
 * @param h
 *   The bandwidth: the percentage of the range to use for clustering.
 */
MeanShiftClusters * meanshift2D(double* s, double * t, int n, double h) {
  MeanShiftClusters * msc = malloc(sizeof(MeanShiftClusters));

  int iter = 200;
  double thr = 0.0001;

  double *a = (double *) malloc(sizeof(double) * n);
  double *b = (double *) malloc(sizeof(double) * n);

  // First calculate the size of the range for each vector
  int i, j;
  double si[2];
  double a_min = INFINITY;
  double b_min = INFINITY;
  double a_max = -INFINITY;
  double b_max = -INFINITY;

  for (i = 0; i < n; i++) {

    if (a_min > s[i]) {
      a_min = s[i];
    }
    if (a_max < s[i]) {
      a_max = s[i];
    }

    if (b_min > t[i]) {
      b_min = t[i];
    }
    if (b_max < t[i]) {
      b_max = t[i];
    }
  }
  si[0] = a_max - a_min;
  si[1] = b_max - b_min;

  // Scale s and t by si. They become a and b. In the original R code,
  // the option to scale was provided in the function arguments. Here it is
  // applied automatically.
  for (i = 0; i < n; i++) {
    a[i] = s[i] / si[0];
    b[i] = t[i] / si[1];
  }

  double **finals = (double **) malloc (sizeof(double*) * n);
  double **savecluster = (double **) malloc (sizeof(double*) * n);
  int ncluster = 0;
  int * cluster_label = (int *) malloc (sizeof (int) * n);
  double cluster_dist[n];
  MeanShiftRep * msr;
  double min_dist = 0;
  double which_min = 0;

  // Iterate through the elements of a & b (i.e the rows of the 2D matrix)
  double x[2];
  for (i = 0; i < n; i++) {

    // Create the point, x.
    x[0] = a[i];
    x[1] = b[i];

    msr = meanshift_rep(a, b, n, x, h, 1e-8, iter);
    finals[i] = (double *) malloc(sizeof(double) * 2);
    finals[i][0] = msr->final[0];
    finals[i][1] = msr->final[1];

    // Initialize the cluster_dist with zeros but only for the number of
    // clusters that we have. The number of clusters is determined by the code
    // further down, so the value of ncluster can change on each iteration if i.
    for (j = 0; j < ncluster; j++) {
      cluster_dist[j] = 0;
    }

    // If we have a cluster then calculate the distance of this point to the
    // cluster and find the one  with the minimum distance.  This will be the
    // cluster to which the point, x, belongs
    if (ncluster >= 1) {
      min_dist = INFINITY;
      for (j = 0; j < ncluster; j++) {
        double t[2];
        t[0] = savecluster[j][0] - finals[i][0];
        t[1] = savecluster[j][1] - finals[i][1];
        cluster_dist[j] = euclidian_norm(t, 2) / euclidian_norm(savecluster[j], 2);
        if (min_dist > cluster_dist[j]) {
          min_dist = cluster_dist[j];
          which_min = j;
        }
      }
    }
    // If we have no clusters or the minimum distance is greater than the
    // threshold then perform the following.
    if (ncluster == 0 || min_dist > thr) {
      savecluster[ncluster] = (double *) malloc(sizeof(double) * 2);
      savecluster[ncluster][0] = finals[i][0];
      savecluster[ncluster][1] = finals[i][1];
      ncluster = ncluster + 1;
      cluster_label[i] = ncluster;
    }
    // If we have a cluster and the minimum distance is less than the threshold
    // then add a label to the point
    else {
      cluster_label[i] = which_min + 1;
    }
    free_msr(msr);
  }
  // Free allocated memory
  for (i = 0; i < n; i++) {
    free(finals[i]);
  }
  free(finals);

  // Find the nearest cluster center in euclidean distance
  /*
  int closest_label[n];
  for (i = 0; i < n; i++){
    // Create the point, x.
    x[0] = a[i];
    x[1] = b[i];
    int md = minimal_dist(savecluster, ncluster, x);
    closest_label[i] = md;
  }*/

  // Count the number of elements in each cluster
  msc->sizes = (int *) malloc(sizeof (int) * ncluster);
  for (i = 0; i < ncluster; i++) {
    msc->sizes[i] = 0;
  }
  for (i = 0; i < n; i++) {
    int cluster = cluster_label[i] - 1;
    msc->sizes[cluster]++;
  }

  // Organize the points into their corresponding clusters
  int counts[ncluster];
  msc->clusters = (double ***) malloc(sizeof (double **) * ncluster);
  for (i = 0; i < ncluster; i++) {
    msc->clusters[i] = (double **) malloc(sizeof(double *) * msc->sizes[i]);
    counts[i] = 0;
  }
  for (i = 0; i < n; i++) {
    int cluster = cluster_label[i] - 1;
    int index = counts[cluster];
    msc->clusters[cluster][index] = (double *) malloc(sizeof(double) * 2);
    msc->clusters[cluster][index][0] = s[i];
    msc->clusters[cluster][index][1] = t[i];
    counts[cluster]++;
  }

  // return the clusters
  msc->num_clusters = ncluster;
  msc->cluster_label = cluster_label;
  msc->centers = savecluster;
  msc->a = a;
  msc->b = b;
  msc->scaled_by[0] = si[0];
  msc->scaled_by[1] = si[1];
  return msc;
}

/**
 * Frees the memory from a MeanShiftClusters object.
 *
 * @param MeanShiftClusters * msc
 *   An instantiated MeanShiftClusters
 */
void free_msc(MeanShiftClusters * msc) {

  int k, l;

  for(k = 0; k < msc->num_clusters; k++) {
    for (l = 0; l < msc->sizes[k]; l++) {
      free(msc->clusters[k][l]);
    }
    free(msc->clusters[k]);
    free(msc->centers[k]);
  }
  free(msc->sizes);
  free(msc->clusters);
  free(msc->cluster_label);
  free(msc->centers);
  free(msc->a);
  free(msc->b);
  free(msc);
}
/**
 * Frees the memory from a MeanShiftRep object.
 *
 * @param MeanShiftRep * msr
 *   An instantiated MeanShiftRep
 */
void free_msr(MeanShiftRep * msr) {
  int j;
  free(msr->start);
  free(msr->thresh);

  for (j = 0; j < msr->iterations; j++) {
    free(msr->points[j]);
  }
  free(msr->points);
  free(msr->final);
  free(msr);
}
/**
 * Mean shift iterative function (until convergence ...)
 */
MeanShiftRep * meanshift_rep(double* a, double * b, int n, double * x, double h, double thresh, int iter) {
  int d = 2;
  int s = -1;
  int j = 0;
  double *x0 = malloc(sizeof(double) * 2);
  double *xt = malloc(sizeof(double) * 2);
  double **M;
  double *th = malloc(sizeof(double) * iter);
  double *m;
  double mx[d];
  MeanShiftRep * msr = malloc(sizeof(MeanShiftRep));

  // Copy the original x into x0 which keeps the start point
  x0[0] = x[0];
  x0[1] = x[1];

  // Copy the point into the xt for iteration.
  xt[0] = x[0];
  xt[1] = x[1];

  M = (double **) malloc(sizeof(double *) * iter);
  for (j = 0; j < iter; j++) {
    m = meanshift_base(a, b, n, xt, h);
    M[j] = m;
    mx[0] = m[0] - xt[0];
    mx[1] = m[1] - xt[1];
    th[j] = euclidian_norm(mx, 2) / euclidian_norm(xt, 2);

    // on the first iteration we free the xt variable
    if (j == 0) {
      free(xt);
    }

    if (th[j] < thresh){
      s = j;
      break;
    }

    // the m becomes the new x
    xt = m;
  }
  if (s == -1) {
    s = iter - 1;
  }

  msr->points = M;
  msr->iterations = s + 1;
  msr->start = x0;
  msr->thresh = th;
  msr->final = malloc(sizeof(double) * 2);
  msr->final[0] = m[0];
  msr->final[1] = m[1];

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
double minimal_dist(double **x, int n, double *y) {
  int d = 2;
  int i;
  double * dv = distance_vector(x, n, y);
  double s[n];
  double min_s = INFINITY;

  //int min_i = 0;

  for (i = 0; i < n; i++) {
    s[i] = sqrt(d) * dv[i];
    if (s[i] < min_s) {
      min_s = s[i];
      //min_i = i;
    }
  }
  free(dv);

  // Return the minimum distance and the position.
//  double *out = (double *) malloc(sizeof(double) * 2);
//  out[0] = min_s;
//  out[1] = min_i;
  return min_s;
}

/**
 * Computes all distances between a set of vectors and another vector.
 * (up to the constant sqrt(d))
 */
double * distance_vector(double **x, int n, double *y) {
  int i;
  int d = 2;

  // The original R code looks like the following.
  // out <- sqrt(as.vector(rowMeans(X^2) + mean(y^2) - 2 * X %*% y/n)))
  // the %*% means matrix multiplication

  // Calculate the squared mean of each row in our ab Matrix, and the squared
  // mean of y.
  double row_means[n];
  double y_mean = 0;
  for (i = 0; i < n; i++) {
    row_means[i] = (pow(x[i][0], 2) + pow(x[i][1], 2)) / 2.0;
  }
  y_mean = (pow(y[0], 2) + pow(y[1],2)) / d;

  // Get a vector by  performing matrix multiplication between ab and y.
  // Multiply the value by 2 and divide by n.
  double Xy[n];
  for (i = 0; i < n; i++) {
    Xy[i] = (x[i][0] * y[0] + x[i][1] * y[1]);
  }

  // Finish up the distance vector.
  double * out = (double *) malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    out[i] = sqrt(row_means[i] + y_mean - 2 * Xy[i] / d);
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
double * meanshift_base(double *a, double *b, int n, double *x, double h) {
  int i;
  double * g;
  double * ms = (double *) malloc(sizeof(double) * 2);
  double sum_g = 0;
  double sum_ag = 0;
  double sum_bg = 0;
  memset(ms, 0, 2);

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
    double p = pow((x-xi[j]) / h, 2);
    double e = exp(-0.5 * p);
    k1[j] = 0.5 * e;
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

  // Get the 1D profiles for both a and b. We use the corresponding value
  // of x to pass in to the profile1d function (a first, b second).
  pa = profile1d(a, n, x[0], h);
  pb = profile1d(b, n, x[1], h);

  // Calculate the multidimensional profile
  for (i = 0; i < n; i++) {
    k[i] = pa[i] * pb[i];
  }

  free(pa);
  free(pb);

  return k;
}

/**
 * Mean shift clustering bandwidth selection.
 *
 * @param double *s
 * @param double *t
 * @param int n
 *   The size of the s and t arrays.
 * @param float taumin
 * @param fload taumax
 * @param int gridsize
 */
double * meanshift_coverage2D(double *s, double *t, int n) {
  int i, j;

  // Set some default parameters. Perhaps these should be passed in?
  double taumin = 0.02;
  double taumax = 0.5;
  int gridsize = 25;

  double h0 = taumin;
  double h1 = taumax;

  // Create a sequence between the taumin and taumax that are equally spaced.
  double * h = (double *) malloc(sizeof(double) * gridsize);
  double step = (h1 - h0) / (gridsize - 1);

  h[0] = h0;
  for (i = 1; i < gridsize - 1; i++) {
    h[i] = h[i-1] + step;
  }
  h[i] = h1;

  // Initialize the cover matrix
  double cover[2][gridsize];
  for (i = 0; i < 2; i++) {
    for (j = 0; j < gridsize; j++) {
      cover[i][j] = 0;
    }
  }

  for (i = 0; i < gridsize; i++) {
    double new_h0 = h[i];
    // fit <- ms(X, new.h0,  thr = thr, scaled = scaled, plotms = 0,  or.labels=or.labels)
    MeanShiftClusters * fit = meanshift2D(s, t, n, new_h0);

    // Find all clusters with more than 2 members
    // find <- as.numeric(which(table(fit$cluster.label)>2))  # changed 23/05/11
    int find[fit->num_clusters];
    int num_found = 0;
    for (j = 0; j < fit->num_clusters; j++) {
      if (fit->sizes[j] > 2) {
        find[j] = 1;
        num_found++;
      }
      else {
        find[j] = 0;
      }
    }

    // Add the center points for clusters with more than 2 elements into an array.
    double ** Pm = (double **) malloc(sizeof(double *) * num_found);
    int k = 0;
    for (j = 0; j < fit->num_clusters; j++) {
      if (find[j] == 1) {
        Pm[k] = (double *) malloc(sizeof(double) * 2);
        Pm[k][0] = fit->centers[j][0];
        Pm[k][1] = fit->centers[j][1];
        k++;
      }
    }

    cover[0][i] = new_h0;
    cover[1][i] = coverage_raw(fit->a, fit->b, n, Pm, k, new_h0);
    for (j = 0; j < k; j++) {
      free(Pm[j]);
    }
    free(Pm);
    free_msc(fit);
  }

  free(h);
  return select_coverage((double **) &cover, gridsize, 1.0/3.0);
}

/**
 * @param double ** cover
 *   A 2 dimensional vector of doubles. With the first column being a
 *   bandwidth and the second being the coverage.
 * @param int n
 *   The size of the cover array.
 */
double * select_coverage(double ** cover, int n, double smin) {

  int diff1[n], diff2[n];
  int i;

  // Initialize the diff1 and diff2 arrays.
  for (i = 0; i < n; i++) {
    diff1[i] = 0;
    diff2[i] = 0;
  }

  for (i = 1; i < n; i++) {
    diff1[i] = cover[i][1] - cover[i-1][1];
  }
  for (i = 1; i < n - 1; i++) {
    diff2[i] = diff1[i + 1] - diff1[i];
  }

  // select <- select.coverage <- select.2diff <- NULL
  double select[n], select_coverage[n], select_2diff[n];
  int k = 0;
  for (i = 2; i < n - 1; i++) {
    // find the maximum coverage between 0 and i - 1;
    double max = smin;
    int j;
    for (j = 0; j < i - 1; j++) {
      if (max < cover[j][1]) {
        max = cover[j][1];
      }
    }
    if (diff2[i] < 0 && cover[i][1] > max) {
      select[k] = cover[i][0];
      select_coverage[k] = cover[i][1];
      select_2diff[k] = diff2[i];
      k++;
    }
  }

  // reorder the select and select_coverage arraays based on the re-ordering
  // of the select_2diff array
  double * selected = (double *) malloc(sizeof(double) * k);
  double * covered = (double *) malloc(sizeof(double) * k);
  int * order = quickSortOrder((double *)&select_2diff, k);
  for (i = 0; i < k; i++) {
    selected[i] = select[order[k]];
    covered[i] = select_coverage[order[k]];
  }

  // The R code returned covered as part of the select.self.coverage
  // function, but we don't need it just for bandwidth selection.
  free(covered);
  return selected;
}

/**
 * @param double *a
 *   The first scaled data vector returned from meanshift2D
 * @param double *b
 *   The second scaled data vector returned from meanshift2D
 * @param int n
 *   The size of vectors a and b.
 */
double coverage_raw(double * a, double *b, int n, double ** centers, int nc, double tau) {
   double min_distance[n];
   int i;
   for (i = 0; i < n; i++) {
     min_distance[i] = 0;
   }

   for (i = 0; i < n; i++) {
     double * point = (double *) malloc(sizeof(double) * 2);
     point[0] = a[i];
     point[1] = b[i];
     min_distance[i] = minimal_dist(centers, nc, point);
   }
   // Calculate the mean. R code does a weighted mean but we don't need
   // to use a weight.
   double ci = 0;
   int num_less = 0;
   for (i = 0; i < n; i++) {
     if (min_distance[i] <= tau) {
       ci += 1;
     }
   }
   ci = ci / n;
   return ci;
 }
