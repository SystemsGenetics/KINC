#include "royston.h"

/**
 * Calculates the Royston H test for bivariate normality.
 *
 * @param double * a
 * @param double * b
 * @param int n
 *   The size of a and b.
 */
double royston2D(double* a, double * b, int n) {

  // The cols variable is the number of genes. Because this is bivariate it
  // will always be 2.  The rows is the number of sample measurements per gene.
  int cols = 2;
  int rows = n;

  // Square and cube the number of rows. Used for later calculations
  int rows2 = pow(rows, 2);
  int rows3 = pow(rows, 3);

  // Create a new vector z
  // z <- matrix(nrow=p, ncol = 1)
  double z[cols];
  double w[cols];

  // We must have at least 3 rows.
  if (rows <= 3) {
    // TODO: should this raise an error?
    return NAN;
  }
  if (rows > 2000) {
    // stop("n must be less than 2000")
    return NAN;
  }
  // If we have between four and 11 rows
  else if (rows >= 4 && rows <= 11) {
    double g = -2.273 + 0.459 * rows;
    double m = 0.5440 - 0.39978 * rows + 0.025054 * rows2 - 0.0006714 * rows3;
    double s = exp(1.3822 - 0.77857 * rows + 0.062767 * rows2 - 0.0020322 * rows3);

    // Get the measure of kurtosis.
    double ak = kurtosis(a, rows);
    double bk = kurtosis(b, rows);

    // Shapiro-Francia test is better for leptokurtic samples.
    if (ak > 3) {
      w[0] = sfrancia(a, rows);
    }
    // Shapiro-Wilk test is better for platykurtic samples.
    else {
      int ifault = 0;
      double W = 0, pw;
      swilk(a, rows, &W, &pw, &ifault);
      w[0] = W;
      if (ifault > 0 && ifault != 7) {
        // error("ifault=%d. This should not happen", ifault);
        // TODO: deal with this error
        return NAN;
      }
    }
    // Shapiro-Francia test is better for leptokurtic samples.
    if (bk > 3) {
      w[1] = sfrancia(b, rows);
    }
    // Shapiro-Wilk test is better for platykurtic samples.
    else {
      int ifault = 0;
      double W = 0, pw;
      swilk(b, rows, &W, &pw, &ifault);
      w[1] = W;
      if (ifault > 0 && ifault != 7) {
        // error("ifault=%d. This should not happen", ifault);
        // TODO: deal with this error
        return NAN;
      }
    }

    z[0] = (-log(g - (log(1 - w[0]))) - m) / s;
    z[1] = (-log(g - (log(1 - w[1]))) - m) / s;
  }
  // If we have between 12 and 2000 rows
  else if (rows >= 12 && rows <= 2000) {
    double x = log(rows);
    double g = 0;
    double m = -1.5861 - 0.31082 * x - 0.083751 * rows2 + 0.0038915 * rows3;
    double s = exp(-0.4803 -0.082676 * x + 0.0030302 * rows2);

    // Get the measure of kurtosis.
    double ak = kurtosis(a, rows);
    double bk = kurtosis(b, rows);

    // Shapiro-Francia test is better for leptokurtic samples.
    if (ak > 3) {
      w[0] = sfrancia(a, rows);
    }
    // Shapiro-Wilk test is better for platykurtic samples.
    else {
      int ifault = 0;
      double W = 0, pw;
      swilk(a, rows, &W, &pw, &ifault);
      w[0] = W;
      if (ifault > 0 && ifault != 7) {
        //error("ifault=%d. This should not happen", ifault);
        // TODO: deal with this error
        return NAN;
      }
    }
    // Shapiro-Francia test is better for leptokurtic samples.
    if (bk > 3) {
      w[1] = sfrancia(b, rows);
    }
    // Shapiro-Wilk test is better for platykurtic samples.
    else {
      int ifault = 0;
      double W = 0, pw;
      swilk(b, cols, &W, &pw, &ifault);
      w[1] = W;
      if (ifault > 0 && ifault != 7) {
        //error("ifault=%d. This should not happen", ifault);
        // TODO: deal with this error
        return NAN;
      }
    }

    z[0] =  ((log(1 - w[0])) + g - m) / s;
    z[1] =  ((log(1 - w[1])) + g - m) / s;
  }
  else {
    // stop("n is not in the proper range")
    return NAN;
  }

//  double u = 0.715;
//  double v = 0.21364 + 0.015124 * pow(log(rows), 2) - 0.0018034 * pow(log(rows), 3);
//  double l = 5;

  // Get the correlation of a and b
  double pcc = gsl_stats_correlation(a, 1, b, 1, rows);

  // Equivalent degrees of freedom
  double edf = cols / (1 + (cols - 1) * pcc);

  double res[cols];
  int j;
  for (j = 0; j < cols; j++) {
    double pn = pnorm(-z[j], 0, 1, TRUE, FALSE);
    double qn = qnorm(pn / 2, 0, 1, TRUE, FALSE);
    res[j] = pow(qn, 2);
  }

  double RH = (edf * (res[0] + res[1])) / cols;
  //double pv = pchisq(RH, edf, lower.tail = FALSE);
  double pv = gsl_cdf_chisq_P(RH, edf);

  return pv;
}
