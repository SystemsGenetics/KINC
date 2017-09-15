#include "kurtosis.h"

/**
 * Computes the estimator of Pearson’s measure of kurtosis.
 *
 * Kurtosis is a measure of the "peakedness" of the probability distribution.
 * The code below was adapted from the kurtosis function of the moments R
 * package: http://cran.r-project.org/web/packages/moments/
 *
 * @param float *x
 *   A numerical vector.
 * @param int n
 *   The size of the vector
 */
float kurtosis(float *x, int n) {
  int i;
  float sum = 0;
  float mean_x = 0;

  // Calculate mean
  for (i = 0; i < n; i++) {
    sum += x[i];
  }
  mean_x = sum / n;

  // Calculate sum((x-mean(x))^4) and sum((x-mean(x))^2)
  float r = 0;
  float q = 0;
  for (i = 0; i < n; i++) {
    r += pow(x[i] - mean_x,  4);
    q += pow(x[i] - mean_x,  2);
  }

  float k = n *r / pow(q, 2);
  return k;
}
