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

  // Create the clean arrays.
  this->clean();
}
/**
 *
 */
PairWiseSet::PairWiseSet(double *a, double *b, int n, int i, int j) {
  this->gene1 = i;
  this->gene2 = j;

  this->n_orig = n;
  this->x_orig = a;
  this->y_orig = b;

  this->x_clean = NULL;
  this->y_clean = NULL;
  this->n_clean = NAN;
  this->samples = NULL;

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
 */
void PairWiseSet::clean() {
  // Create the vectors that will contain the sample measurements that are
  // not empty (e.g. NAN or INF).
  this->x_clean = (double *) malloc(sizeof(double) * this->n_orig);
  this->y_clean = (double *) malloc(sizeof(double) * this->n_orig);

  // Create the samples array.
  this->samples = (int *) malloc(sizeof(int) * this->n_orig);

  // Iterate through both data arrays. If either of them have an NAN or INF
  // then remove that sample from both array.
  int n = 0;
  for (int i = 0; i < this->n_orig; i++) {
    if (isnan(this->x_orig[i]) || isnan(this->y_orig[i]) || isinf(this->x_orig[i]) || isinf(this->y_orig[i])) {
      // This sample cannot be used, so set the sample to 0 and skip.
      this->samples[i] = 0;
      continue;
    }
    // This sample can be used so set the sample to one and keep the values.
    this->x_clean[n] = this->x_orig[i];
    this->y_clean[n] = this->y_orig[i];
    this->samples[i] = 1;
    n++;
  }
  // Set the final size of the cleaned arrays.
  this->n_clean = n;
}
