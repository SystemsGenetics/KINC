#include "sfrancia.h"
/**
 * Shapiro-Francia normality test
 */
void sfrancia(float *vector, int n, float *w, float *pw, int *ifault) {

  if ((n < 5 || n > 5000)) {
    *ifault = 1;
    char message[100] = "You must have between 3 and 500 samples for Shapiro Francia normality test.";
    handle_warning(message);
    return;
  }

  // Remove missing values from x

  // Create a copy of the vector and sort it.
  float * x = (float *) malloc(sizeof(float) * n);
  memcpy(x, vector, sizeof(float) * n);

  // Sort the incoming vector
  quickSortF(x, n);

  // Get equally spaced points between 0 and 1 and then determine
  // their probabilities from a normal distribution
  float a = 3.0/8.0;
  float * norm_points = ppoints(n, a);

  // Get the value of the probabilities using x
  float y[n];
  int i = 0;
  for (i = 0; i < n; i++) {
    y[i] = qnorm(norm_points[i], 0, 1, FALSE, FALSE);
  }
  free(norm_points);

  // float W = cor(x, y)^2;
  float pcc = gsl_stats_float_correlation(x, 1, y, 1, n);
  float W = pow(pcc, 2);
  float u = log(n);
  float v = log(u);
  float mu = -1.2725 + 1.0521 * (v - u);
  float sig = 1.0308 - 0.26758 * (v + 2/u);
  float z = (log(1.0 - W) - mu) / sig;
  float pval = pnorm(z, 0, 1, FALSE, FALSE);

  // Return the probabiliy
  *w = W;
  *pw = pval;

  free(x);
}
