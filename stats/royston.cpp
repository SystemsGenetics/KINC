#include "royston.h"

/**
 * Calculates the Royston H test for bivariate normality.
 *
 * This function was adapted from the R function royston.test which is
 * available in the royston package at
 * http://cran.r-project.org/web/packages/royston/index.html.
 *
 * Unlike the R code, this function has been adapted to only deal with
 * two vectors rather than a matrix.
 *
 * @param float * a
 * @param float * b
 * @param int n
 *   The size of a and b.
 * @param float * pcc
 *   The pearson's correlation value. It is set during function run
 */
float roystonH(float* a, float * b, int n, float *pcc) {

  // The cols variable is the number of genes. Because this is bivariate it
  // will always be 2.  The rows is the number of sample measurements per gene.
  int cols = 2;
  int rows = n;

  // Create a new vector z
  float z[cols];

  // Initialize the pcc value
  *pcc = NAN;

  // We must have at least 3 rows.
  if (rows <= 3) {
    char message[100] = "You must have at least 3 samples for Royston's H test.";
    handle_warning(message);
    return NAN;
  }
  if (rows > 2000) {
    char message[100] = "You must have no more than 2000 samples for Royston's H test.";
    handle_warning(message);
    return NAN;
  }
  // If we have between four and 11 rows
  else if (rows >= 4 && rows <= 11) {
    float x = n;
    float g = -2.273 + 0.459 * x;
    float m = 0.5440 - 0.39978 * x + 0.025054 * pow(x, 2) - 0.0006714 * pow(x, 3);
    float s = exp(1.3822 - 0.77857 * x + 0.062767 * pow(x, 2) - 0.0020322 * pow(x, 3));

    int i = 0;
    for (i = 0; i < cols; i++) {
      float k = 0;
      int ifault = 0;
      float W = 0, pw;
      float *select;

      if (i == 0) {
        select = a;
      }
      else {
        select = b;
      }

      k = kurtosis(select, rows);
      if (k > 3) {
        sfrancia(select, rows, &W, &pw, &ifault);
      }
      else {
        swilk(select, rows, &W, &pw, &ifault);
      }
      if (ifault > 0 && ifault != 7) {
        char message[100] = "Normality test failed while calculating Royston's H test.";
        handle_warning(message);
        return NAN;
      }
      z[i] = (-log(g - (log(1 - W))) - m)/s;
    }
  }
  // If we have between 12 and 2000 rows
  else if (rows >= 12 && rows <= 2000) {
    float x = log(rows);
    float g = 0;
    float m = -1.5861 - 0.31082 * x - 0.083751 * pow(x,2) + 0.0038915 * pow(x,3);
    float s = exp(-0.4803 -0.082676 * x + 0.0030302 * pow(x,2));

    int i = 0;
    for (i = 0; i < cols; i++) {
      float k = 0;
      int ifault = 0;
      float W = 0, pw;
      float *select;

      if (i == 0) {
        select = a;
      }
      else {
        select = b;
      }

      k = kurtosis(select, rows);
      if (k > 3) {
        sfrancia(select, rows, &W, &pw, &ifault);
      }
      else {
        swilk(select, rows, &W, &pw, &ifault);
      }
      if (ifault > 0 && ifault != 7) {
        char message[100] = "Normality test failed while calculating Royston's H test.";
        handle_warning(message);
        return NAN;
      }
      z[i] = ((log(1 - W)) + g - m) / s;
    }
  }
  else {
    // stop("n is not in the proper range")
    return NAN;
  }

  float u = 0.715;
  float v = 0.21364 + 0.015124 * pow(log(rows), 2) - 0.0018034 * pow(log(rows), 3);
  float l = 5;

  // Get the correlation of a and b
  *pcc = gsl_stats_float_correlation(a, 1, b, 1, rows);

  // Transformed PCC value
  float NC = pow(*pcc, l) * (1.0 - (u * pow(1.0 - *pcc, u)) / v);

  // Calculate the % Total. In the R code this was
  //   T = sum(sum(NC)) - p
  // but, because we only have two vectors it was adjusted as belo.
  float T = 2 + 2 * NC - cols;

  // Calculate the average correlation
  float mC = T / (pow(cols, 2) - cols);

  // Equivalent degrees of freedom
  float edf = cols / (1.0 + (cols - 1.0) * mC);

  float res[cols];
  int j;
  float sum_res = 0;
  for (j = 0; j < cols; j++) {
    float pn = pnorm(-z[j], 0, 1, TRUE, FALSE);
    float qn = qnorm(pn / 2, 0, 1, TRUE, FALSE);
    res[j] = pow(qn, 2);
    sum_res += res[j];
  }

  float RH = (edf * (sum_res)) / cols;
  //float pv = pchisq(RH, edf, lower.tail = FALSE);
  float pv = 1 - gsl_cdf_chisq_P(RH, edf);

  return pv;
}
