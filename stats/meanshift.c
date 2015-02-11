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

  free(a);
  free(b);

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
void meanshift_coverage2D(double *s, double *t, int n) {
  int i, j;

  // Set some default parameters. Perhaps these should be passed in?
  float taumin = 0.02;
  float taumax = 0.5;
  int gridsize = 25;

  float Pm = NULL;
  float h0 = taumin;
  float h1 = taumax;

  // Create a sequence between the taumin and taumax that are equally spaced.
  float * h = (float *) malloc(sizeof(float) * gridsize);
  int step = (h1 - h0) / (gridsize - 2);

  h[i] = h0;
  for (i = 1; i < gridsize - 1; i++) {
    h[i] = h[i-1] + step;
  }
  h[i] = h1;

  // Initialize the cover matrix
  float cover[2][gridsize];
  for (i = 0; i < 2; i++) {
    for (j = 0; j < gridsize; j++) {
      cover[i][j] = 0;
    }
  }

  for (i = 0; i < gridsize; i++) {
    double new_h0 = h[i];
    // fit <- ms(X, new.h0,  thr = thr, scaled = scaled, plotms = 0,  or.labels=or.labels)
    MeanShiftClusters * msc = meanshift2D(s, t, n, new_h0);

    // Find all clusters with more than 2 members
    // find <- as.numeric(which(table(fit$cluster.label)>2))  # changed 23/05/11
    int find[msc->num_clusters];
    for (j = 0; j < msc->num_clusters; j++) {
      if (msc->sizes[j] > 2) {
        find[j] = 1;
      }
      else {
        find[j] = 0;
      }
    }

    Pm[[i]] <- fit$cluster.center[find,]
    if (!is.matrix(Pm[[i]])){
      Pm[[i]]<- matrix(Pm[[i]],nrow=1, dimnames=list(dimnames(fit$cluster.center)[[1]][find] ,NULL))
    }
    if (!cluster) {
        cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]],  new.h0, plot.type = 0, print=print)[1:2])
    } else {
        cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]], new.h0, plot.type = 0, label = fit$cluster.label, print=print)[1:2])
    }

  }
  select <- select.self.coverage(self = cover,
      smin = 1/3, plot.type = plot.type)
  result <- list(self.coverage.curve = cover, select = select,
      type = "ms")
  class(result) <- "self"
  result
}

coverage.raw <-function(X, vec, tau, weights=1, plot.type="p", print=FALSE, label=NULL,...){
     X<- as.matrix(X)
     p <- dim(vec)[1]
     n <- dim(X)[1]
     d <- dim(X)[2]

     if (n %% length(weights) !=0){
         stop("The length of the vector of weights is not a multiple of the sample size.")
     } else {
      weights <-rep(weights, n %/%   length(weights))
     }

     min.distance  <- rep(0,n)
     ins.distance <- rep(1,n)

     for (i in 1: n){
        #if (i %%10==0){ print(i)}
        if (is.null(label)){
             min.distance[i] <- mindist(vec, X[i,])$mindist
        } else {
            if (!is.matrix(vec)){vec<-matrix(vec, nrow=1)}
            if (as.character(label[i]) %in%  dimnames(vec)[[1]]){
                  min.distance[i] <- mindist(matrix(vec[as.character(label[i]),],nrow=1), X[i,])$mindist   # 11/10/10 experimental, for clustering
            } else {
                min.distance[i]<- tau+1   #01/11/10 if data point not allocatable to any centre.
            }
        }
        ins.distance[i] <- (min.distance[i] <= tau) #indicator for being inside/outside the tube
           }
     ci<- weighted.mean(min.distance <= tau, w=weights)
     if (plot.type %in% c("p","l")) {
        plot(X, col=ins.distance+1,...)
        if (plot.type=="p"){points(vec, col=3,lwd=2 )} else if (plot.type=="l"){lines(vec, col=3,lwd=2 )}
        }
     if (print){print(c(tau,ci))}
     return(list(tau=tau, coverage= ci, min=min.distance, inside= ins.distance))
 }

select.self.coverage <-
function (self,  smin, plot.type = "o", plot.segments=NULL)
{
    if (class(self) == "self") {
        cover <- self$self.coverage.curve
    }
    else {
        cover <- self
    }
    if (missing(smin)) {
        if (class(self) == "self") {
            smin <- switch(self$type, lpc = 2/3, ms = 1/3)
        }
        else stop("Please specify `smin' argument.")
    }
    n <- dim(cover)[1]
    diff1 <- diff2 <- rep(0, n)
    diff1[2:n] <- cover[2:n, 2] - cover[1:(n - 1), 2]
    diff2[2:(n - 1)] <- diff1[3:n] - diff1[2:(n - 1)]
    select <- select.coverage <- select.2diff <- NULL

    if (plot.type != 0) {
        plot(cover, type = plot.type, xlab = "h", ylab = "S(h)",
            ylim = c(0, 1))
    }
    for (i in (3:(n - 1))) {
        if (diff2[i] < 0 && cover[i, 2] > max(smin, cover[1:(i -
            1), 2])) {
            select <- c(select, cover[i, 1])
            select.coverage <- c(select.coverage, cover[i,2])
            select.2diff <- c(select.2diff, diff2[i])
            #if (plot.type != 0) {
            #    segments(cover[i, 1], 0, cover[i, 1], cover[i,
            #      2], col = scol[i], lty = slty[i], lwd=slwd[i])
            #}
        }
    }

    selected <-select[order(select.2diff)]
    covered<- select.coverage[order(select.2diff)]

    if (plot.type != 0) {
      d<-length(selected)
      slty <- slwd <-scol<-rep(0,d)
      if (is.null(plot.segments)){
         slty[1:3]<- c(1,2,3)
         slwd[1:3] <-c(2,1,1)
         scol[1:3] <- c(3,3,3)
      } else {
          r<-max(length(plot.segments$lty), length(plot.segments$lwd), length(plot.segments$col))
            slty[1:r] <- slwd[1:r] <-scol[1:r]<-1
            if (length(plot.segments$lty)>0){ slty[1:length(plot.segments$lty)]<- plot.segments$lty}
             if (length(plot.segments$lwd)>0){slwd[1:length(plot.segments$lwd)]<- plot.segments$lwd}
            if (length(plot.segments$col)>0){ scol[1:length(plot.segments$col)]<- plot.segments$col}
       }
      for (j in 1:d){
           segments(selected[j], 0, selected[j], covered[j], col = scol[j], lty = slty[j], lwd=slwd[j])
         }
      }

    return(list("select"=selected, "select.coverage"=covered, "select.2diff"=select.2diff[order(select.2diff)]))
}
