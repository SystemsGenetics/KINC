#include "sfrancia.h"
/**
 * Shapiro-Francia normality test
 */
void sfrancia(double *vector, int n, double *w, double *pw, int *ifault) {

  if ((n < 5 || n > 5000)) {
    *ifault = 1;
    char message[100] = "You must have between 3 and 500 samples for Shapiro Francia normality test.";
    handle_warning(message);
    return;
  }

  // Remove missing values from x

  // Create a copy of the vector and sort it.
  double * x = (double *) malloc(sizeof(double) * n);
  memcpy(x, vector, sizeof(double) * n);

  // Sort the incoming vector
  quickSortD(x, n);

  // Get equally spaced points between 0 and 1 and then determine
  // their probabilities from a normal distribution
  double a = 3.0/8.0;
  double * norm_points = ppoints(n, a);

  // Get the value of the probabilities using x
  double y[n];
  int i = 0;
  for (i = 0; i < n; i++) {
    y[i] = qnorm(norm_points[i], 0, 1, FALSE, FALSE);
  }
  free(norm_points);

  // double W = cor(x, y)^2;
  double pcc = gsl_stats_correlation(x, 1, y, 1, n);
  double W = pow(pcc, 2);
  double u = log(n);
  double v = log(u);
  double mu = -1.2725 + 1.0521 * (v - u);
  double sig = 1.0308 - 0.26758 * (v + 2/u);
  double z = (log(1.0 - W) - mu) / sig;
  double pval = pnorm(z, 0, 1, FALSE, FALSE);

  // Return the probabiliy
  *w = W;
  *pw = pval;

  free(x);
}
