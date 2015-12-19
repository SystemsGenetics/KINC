#include "PairWiseSet.h"

/**
 *
 */
PairWiseSet::PairWiseSet(EMatrix * ematrix, int i, int j) {
  this->gene1 = i;
  this->gene2 = j;

  this->n_orig = ematrix->getNumSamples();
  this->x_orig = ematrix->getRow(this->gene1);
  this->y_orig = ematrix->getRow(this->gene2);

  this->x_clean = NULL;
  this->y_clean = NULL;
  this->n_clean = NAN;
  this->samples = NULL;
  this->threshold= NAN;

  // Create the clean arrays.
  this->clean();
}
/**
 *
 */
PairWiseSet::PairWiseSet(EMatrix * ematrix, int i, int j, double th) {
  this->gene1 = i;
  this->gene2 = j;

  this->n_orig = ematrix->getNumSamples();
  this->x_orig = ematrix->getRow(this->gene1);
  this->y_orig = ematrix->getRow(this->gene2);

  this->x_clean = NULL;
  this->y_clean = NULL;
  this->n_clean = NAN;
  this->samples = NULL;
  this->threshold= th;

  // Create the clean arrays.
  this->clean();
}
/**
 *
 */
PairWiseSet::PairWiseSet(double *a, double *b, int n, int i, int j) {
  PairWiseSet(a, b, n, i, j, NAN);
}
/**
 *
 */
PairWiseSet::PairWiseSet(double *a, double *b, int n, int i, int j, double th) {
  this->gene1 = i;
  this->gene2 = j;

  this->n_orig = n;
  this->x_orig = a;
  this->y_orig = b;

  this->x_clean = NULL;
  this->y_clean = NULL;
  this->n_clean = NAN;
  this->samples = NULL;
  this->threshold= th;

  // Create the clean arrays.
  this->clean();
}
/**
 *
 */
PairWiseSet::~PairWiseSet(){
  if (this->samples) {
    free(this->samples);
  }
  free(x_clean);
  free(y_clean);
}
/**
 * Removes the NA's from the sample and set the samples array.
 *
 * @param th
 *   The threshold for marking values as NA.  Any value less than or
 *   equal to this value.
 */
void PairWiseSet::clean() {
  // Create the vectors that will contain the sample measurements that are
  // not empty (e.g. NAN or INF).
  this->x_clean = (double *) malloc(sizeof(double) * this->n_orig);
  this->y_clean = (double *) malloc(sizeof(double) * this->n_orig);

  // Create the samples array.
  this->samples = (int *) malloc(sizeof(int) * this->n_orig);


  // Make the positions in the original vectors that have a
  // value less than or equal to the threshold as NA. We'll also
  // keep track of these positions so we can reclean the pariwise set
  // and mark these as being removed due to the threshold.
  int th_index[n_orig];
  for (int i = 0; i < n_orig; i++) {
    th_index[i] = 0;
  }
  if (!isnan(threshold)) {
    for (int i = 0; i < n_orig; i++) {
      if (x_orig[i] <= threshold) {
        th_index[i] = 1;
        x_orig[i] = NAN;
      }
      if (y_orig[i] <= threshold) {
        th_index[i] = 1;
        y_orig[i] = NAN;
      }
    }
  }

  // Iterate through both data arrays. If either of them have an NAN or INF
  // then remove that sample from both array.
  int n = 0;
  for (int i = 0; i < this->n_orig; i++) {
    // If this sample was removed because one of the x or y values was less than
    // the threshold then set the sample to 6.
    if (th_index[i] == 1) {
      this->samples[i] = 6;
      continue;
    }
    // If the sample is NAN then mark it as missing with a 9.
    if (isnan(this->x_orig[i]) || isnan(this->y_orig[i]) || isinf(this->x_orig[i]) || isinf(this->y_orig[i])) {
      // This sample cannot be used, so set the sample to 0 and skip.
      this->samples[i] = 9;
      continue;
    }
    // Otherwise, this sample can be used so set the sample to 1 and copy
    // the value into the clean vectors.
    this->x_clean[n] = this->x_orig[i];
    this->y_clean[n] = this->y_orig[i];
    this->samples[i] = 1;
    n++;
  }
  // Set the final size of the cleaned arrays.
  this->n_clean = n;
}

/**
 *  Masks outliers
 */
void PairWiseSet::maskOutliers() {

  // Create and initialize vectors for testing outliers.
  double cx[n_clean];
  double cy[n_clean];
  for (int j = 0; j < n_clean; j++) {
    cx[j] = x_clean[j];
    cy[j] = y_clean[j];
  }

  // Discover any outliers
  Outliers * outliersCx = NULL;
  Outliers * outliersCy = NULL;
  outliersCx = outliers_iqr(cx, n_clean, 1.5);
  outliersCy = outliers_iqr(cy, n_clean, 1.5);

  // Mask any outliers in the samples array.
  for (int i = 0; i < n_orig; i++) {
    for (int j = 0; j < outliersCx->n; j++) {
      if (x_orig[i] == outliersCx->outliers[j]) {
        samples[i] = 7;
      }
    }
    for (int j = 0; j < outliersCy->n; j++) {
      if (y_orig[i] == outliersCy->outliers[j]) {
        samples[i] = 7;
      }
    }
  }
  if (outliersCx) {
    free(outliersCx->outliers);
    free(outliersCx);
  }
  if (outliersCy) {
    free(outliersCy->outliers);
    free(outliersCy);
  }

  // Now recreate the clean arrays but with the outliers missing.
  double * nx_clean = (double *) malloc(sizeof(double) * this->n_orig);
  double * ny_clean = (double *) malloc(sizeof(double) * this->n_orig);
  int n = 0;
  for (int i = 0; i < n_orig; i++) {
    if (samples[i] == 1) {
      nx_clean[n] = x_orig[i];
      ny_clean[n] = y_orig[i];
      n++;
    }
  }
  n_clean = n;
  free(x_clean);
  free(y_clean);
  x_clean = nx_clean;
  y_clean = ny_clean;
}
