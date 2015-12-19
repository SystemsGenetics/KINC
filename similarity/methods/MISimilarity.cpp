#include "MISimilarity.h"


MISimilarity::MISimilarity(PairWiseSet * pws, int min_obs, int * samples)
  :PairWiseSimilarity(pws, min_obs, samples) {

  this->mi_bins = 10;
  this->mi_degree = 3;

  strcpy(this->type, "mi");
}

MISimilarity::MISimilarity(PairWiseSet * pws, int min_obs, double mi_bins, double mi_degree)
  :PairWiseSimilarity(pws, min_obs) {

  this->mi_bins = mi_bins;
  this->mi_degree = mi_degree;

  strcpy(this->type, "mi");
}

/**
 * Constructor.
 */
MISimilarity::MISimilarity(PairWiseSet * pws, int min_obs, int * samples, double mi_bins, double mi_degree)
  :PairWiseSimilarity(pws, min_obs, samples) {

  this->mi_bins = mi_bins;
  this->mi_degree = mi_degree;

  strcpy(this->type, "mi");
}

MISimilarity::~MISimilarity() {

}
/**
 * Performs Mutual Information analysis on two arrays.
 *
 * @param double *a
 * @param double *b
 * @param int n
 */
void MISimilarity::run() {
  // Make sure we have the correct number of observations before performing
  // the comparison.
  if (this->n >= this->min_obs) {
    // Calculate the min and max
    double xmin = INFINITY;
    double ymin = INFINITY;
    double xmax = -INFINITY;
    double ymax = -INFINITY;
    for (int i = 0; i < this->n; i++) {
      // calculate the x and y minimum
      if (this->a[i] < xmin) {
        xmin = this->a[i];
      }
      if (this->a[i] > xmax) {
        xmax = this->a[i];
      }
      if (this->b[i] < ymin) {
        ymin = this->b[i];
      }
      if (this->b[i] > ymax) {
        ymax = this->b[i];
      }
    }

    // Make sure that the min and max are not the same.
    if(xmin < xmax && ymin < ymax) {
      score = calculateBSplineMI(this->a, this->b, this->n,
          this->mi_bins, this->mi_degree, xmin, ymin, xmax, ymax);
    }
  }
  else {
    score = NAN;
  }
}

/**
 * Calculates the Mutual Information matrix
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 * @param double ** data
 *   A pointer to a two dimensional array of doubles containing the
 *   expression values
 * @param int * histogram
 *   An integer array created by the init_histogram function and
 *   used for building the histogram.
 *
 */
/*void MISimilarity::calculate_MI(CCMParameters params, double ** data, int * histogram) {
  int j, k, m;             // integers for looping
  char outfilename[1024];  // the output file name
  float max_mi = 0;


  // output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
  int z = (params.rows - 1) / ROWS_PER_OUTPUT_FILE;
  int j_limit;

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./MI", &st) == -1) {
      mkdir("./MI", 0700);
  }

  int total_comps = (params.rows * params.rows) / 2;
  int n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating mutual information...\n");
  for (m = 0; m <= z; m++) {

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./MI/%s.mi%d.bin", params.fileprefix, m);
    printf("Writing file %d of %d: %s... \n", m + 1, z + 1, outfilename);

    FILE * outfile = fopen(outfilename, "wb");

    // calculate the limit on the rows to output based on where we are in the calculation
    if (m < z) {
      j_limit = (m + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      j_limit = params.rows;
    }

    // output the size of the overall matrix
    fwrite(&params.rows, sizeof(params.rows), 1, outfile);

    // determine and output the number of rows that will be stored in the current file
    int i = j_limit - m * ROWS_PER_OUTPUT_FILE;
    fwrite(&i, sizeof(i), 1, outfile);

    for (j = m * ROWS_PER_OUTPUT_FILE; j < j_limit; j++) {
      // lower triangular symmetric matrix (stop when col# > row#)
      for (k = 0; k <= j; k++) {

        n_comps++;
        if (n_comps % 100 == 0) {
          printf("Percent complete: %.2f%%\r", (n_comps/(float)total_comps)*100);
        }

        // build the vectors for calculating MI
        double x[params.cols + 1];
        double y[params.cols + 1];
        double xmin = 9999999;
        double ymin = 9999999;
        double xmax = -9999999;
        double ymax = -9999999;
        int n = 0;
        for (i = 0; i < params.cols; i++) {
          // if either of these elements is missing then don't include the
          // elements from this sample
          if (isnan(data[j][i]) || isnan(data[k][i]) || isinf(data[j][i]) || isinf(data[k][i])) {
            continue;
          }
          // save the x & y vectors
          x[n] = data[j][i];
          y[n] = data[k][i];

          // calculate the x and y minimum
          if (x[n] < xmin) {
            xmin = x[n];
          }
          if (x[n] > xmax) {
            xmax = x[n];
          }
          if (y[n] < ymin) {
            ymin = y[n];
          }
          if (y[n] > ymax) {
            ymax = y[n];
          }
          n++;
        }

        // Calculate Mutual Information if we have enough observations.  We default the
        // MI value to NaN if we do not have the minimum number
        // of observations to do the calculation.  This is better than storing
        // a zero which indicates no correlation.  If we stored zero then
        // we'd be providing a false correlation as no correlation calculation
        // actually occurred.
        double mi = NAN;

        // Make sure we have the right number of observations and that the
        // min and max are not the same.
        if (n >= params.min_obs && xmin < xmax && ymin < ymax) {
          //printf("%d, %d\n", j, k);
          mi = calculateBSplineMI(x, y, n, params.mi_bins, params.mi_degree, xmin, ymin, xmax, ymax);
          //printf("(%d, %d) = %f (%d values). xmin = %f, xmax = %f, ymin = %f, ymax = %f\n", j, k, mi, n, xmin, xmax, ymin, ymax);
        }
        float outmi = (float) mi;
        if (outmi > max_mi) {
          max_mi = outmi;
        }
        fwrite(&outmi, sizeof(float), 1, outfile);


        // if the historgram is turned on then store the value in the correct bin
//        if (mi < 1 && mi > -1) {
//          if (mi < 0) {
//            mi = -mi;
//          }
//          histogram[(int)(mi * HIST_BINS)]++;
//        }
      }
    }
    fclose(outfile);
  }
  printf("\nmax mi: %f\n", max_mi);
}*/

/**
 * @param double *x
 *   A vector of expression values
 * @param double *y
 *   A vector of expression values the same size as x
 * @param int n
 *   The number of data points to fit (or, the size of vectors x and y)
 * @param int m
 *   The number of bins or break points
 * @param int k
 *   The B-splines order or degree of the polynomial.
 * @param int xmin
 * @param int ymin
 * @param int xmax
 * @param int ymax
 *
 */
double MISimilarity::calculateBSplineMI(double *v1, double *v2, int n, int m, int k, double xmin, double ymin, double xmax, double ymax) {

  // used for iterating through loops
  int i, j, q;
  // holds the number of histogram "bins" used to represent the
  // probability distribution for each gene
  int ncoeffs;
  // the final Mutual Information value
  double mi;
  // the observations (e.g. gene x, gene y)
  gsl_vector *x, *y;
  // the spline coefficent vectors
  gsl_vector *Bx, *By;
  // all spline coefficient vecotrs stored in a matrix
  gsl_matrix *BX, *BY;
  // probability distribution vectors
  gsl_vector *px, *py;
  gsl_matrix *pxy;
  // the B-spline workspace structure
  gsl_bspline_workspace * bw;

  // Copy the values into a gsl_vector
  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, i, v1[i]);
    gsl_vector_set(y, i, v2[i]);
  }

  // Allocate a bspline workspace and add the knots.
  bw = gsl_bspline_alloc(k, m);
  // a uniform distribution of knots between zero and 1.  For example,
  // where k = 3 and m =5 the knots will be: 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1
  gsl_bspline_knots_uniform(0, 1, bw);
/*  gsl_vector * knots;
  knots = gsl_vector_alloc(m);
  gsl_vector_set(knots, 0, 0);
  gsl_vector_set(knots, 1, 0);
  gsl_vector_set(knots, 2, 1);
  gsl_vector_set(knots, 3, 2);
  gsl_vector_set(knots, 4, 3);
  gsl_bspline_knots(knots, bw);*/

  // normalize the observations so they fit within the domain of the knot vector [0, 1]
  gsl_vector_add_constant(x, - xmin);
  gsl_vector_scale(x, 1/(xmax - xmin));
  gsl_vector_add_constant(y, - ymin);
  gsl_vector_scale(y, 1/(ymax - ymin));

  // calculate MI
  ncoeffs = gsl_bspline_ncoeffs(bw);
  Bx = gsl_vector_alloc(ncoeffs);
  By = gsl_vector_alloc(ncoeffs);
  BX = gsl_matrix_alloc(n, ncoeffs);
  BY = gsl_matrix_alloc(n, ncoeffs);

  // iterate through each value of x and y and retrieve
  // the co-efficients for the splines these will serve as weights
  // for the probability distribution histograms used by MI
  for (i = 0; i < n; ++i) {
    gsl_bspline_eval(gsl_vector_get(x, i), Bx, bw);
    gsl_bspline_eval(gsl_vector_get(y, i), By, bw);
    // save these co-efficients in a matrix for later use
    for (j = 0; j < ncoeffs; j++) {
      gsl_matrix_set(BX, i, j, gsl_vector_get(Bx, j));
      gsl_matrix_set(BY, i, j, gsl_vector_get(By, j));
    }
  }

  // construct the probability distribution vectors used for MI for each gene
  px = gsl_vector_alloc(ncoeffs);  // probability distribution of x
  py = gsl_vector_alloc(ncoeffs);  // probability distribution of y
  for (j = 0; j < ncoeffs; j++) {
    double px_j = 0;  // the probability of x in bin j
    double py_j = 0;  // the probability of y in bin j
    for (i = 0; i < n; ++i) {
      double wx = gsl_matrix_get(BX, i, j);
      double wy = gsl_matrix_get(BY, i, j);
      px_j += wx;
      py_j += wy;
    }
    px_j /= n;
    py_j /= n;
    gsl_vector_set(px, j, px_j);
    gsl_vector_set(py, j, py_j);
  }

  // construct the joint probability distribution
  pxy = gsl_matrix_alloc(ncoeffs, ncoeffs); // joint probability distribution
  for (j = 0; j < ncoeffs; j++) {
    for (q = 0; q < ncoeffs; ++q) {
      double pxy_jq =  0;
      for (i = 0; i < n; ++i) {
        double wxj = gsl_matrix_get(BX, i, j);
        double wyq = gsl_matrix_get(BY, i, q);
        pxy_jq += (wxj * wyq);
      }
      pxy_jq /= n;
      gsl_matrix_set(pxy, j, q, pxy_jq);
    }
  }

  // calculate the shannon entropy for x, y
  double hx = 0; // shannon entropy for x
  double hy = 0; // shannon entropy for y
  for (j = 0; j < ncoeffs; j++) {
    double px_j = gsl_vector_get(px, j);
    double py_j = gsl_vector_get(py, j);
    if (px_j != 0) {
      hx += px_j * log2(px_j);
    }
    if (py_j != 0) {
      hy += py_j * log2(py_j);
    }
  }
  hx = - hx;
  hy = - hy;

  // calculate the shannon entropy for the joint x,y
  double hxy = 0;
  for (j = 0; j < ncoeffs; j++) {
    for (q = 0; q < ncoeffs; ++q) {
      double pxy_jq = gsl_matrix_get(pxy, j, q);
      if (pxy_jq != 0) {
        hxy += pxy_jq * log2(pxy_jq);
      }
    }
  }
  hxy = - hxy;

  // calculate the mutual information using the formula
  // MI(A,B) = H(A) + H(B) - H(A,B)
  mi = hx + hy - hxy;

  // free up allocated memory
  gsl_bspline_free(bw);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(Bx);
  gsl_vector_free(By);
  gsl_vector_free(px);
  gsl_vector_free(py);
  gsl_matrix_free(BX);
  gsl_matrix_free(BY);
  gsl_matrix_free(pxy);

  return mi;
}
