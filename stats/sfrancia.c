#include "sfrancia.h"
/**
 * Shapiro-Francia normality test
 */
double sfrancia(double *x, int n) {

  if ((n < 5 || n > 5000)) {
    // stop("sample size must be between 5 and 5000")
    return NAN;
  }

  // DNAME <- deparse(substitute(x))

  // Remove missing values from x

  // Sort X
  // double x <- sort(x[complete.cases(x)])

  // Get equally spaced points between 0 and 1 and then determine
  // their probabilities from a normal distribution
  double a = 3/8;
  double * norm_points = ppoints(n, a);

  // Get the value of the probabilities using x
  double y[n];
  int i = 0;
  for (i = 0; i < n; i++) {
    y[i] = qnorm(norm_points[i], 0, 1, FALSE, FALSE);
  }
  // double W = cor(x, y)^2;
  double pcc = gsl_stats_correlation(x, 1, y, 1, n);
  double W = pow(pcc, 2);
  double u = log(n);
  double v = log(u);
  double mu = -1.2725 + 1.0521 * (v - u);
  double sig = 1.0308 - 0.26758 * (v + 2/u);
  double z = (log(1 - W) - mu) / sig;

  // Return the probabiliy
  return pnorm(z, 0, 1, TRUE, FALSE);
}
