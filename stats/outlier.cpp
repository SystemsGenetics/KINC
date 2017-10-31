#include "outlier.h"

/**
 * Finds outliers in an array of points.
 *
 * Identifies the interquantile range for the provided values and finds
 * values that fall outside of that range.  The code for this function
 * has been adapted from the stats::fivenum and boxplot.stats functions of R.
 *
 * @param float * x
 *   A sorted array of doubles.  There must not be any missing values.
 * @param int n
 *   The size of the array 'x'
 * @param float coef
 *   A coefficient to multiply by the IQR. Default should be 1.5.
 *
 * @return Outliers
 *   An Outliers struct.
 */
Outliers * outliers_iqr(float * x, int n, float coef) {

  int i;
  Outliers * outliers = (Outliers *) malloc(sizeof(Outliers));

  // Sort x
  float * sx = (float *) malloc(sizeof(float) * n);
  for (i =0; i < n; i++) {
    sx[i] = x[i];
  }
  quickSortF(sx, n);

  if (coef < 0) {
    char message[100] = "'coef' must not be negative";
    handle_error(message);
  }

  // Initalize the outliers struct
  outliers->outliers = (float *) malloc(sizeof(float) * n);
  outliers->n = 0;

  // Calculate the Interquantile range of the values provided.
  float n4 = ((int)((n + 3)/2.0)) / 2.0;
  float d[5];
  d[0] = 1;
  d[1] = n4;
  d[2] = (n + 1) / 2.0;
  d[3] = n + 1 - n4;
  d[4] = n;

  float stats[5];
  stats[0] = 0.5 * (sx[(int) d[0] - 1] + sx[(int) (d[0] + 0.5) - 1]);
  stats[1] = 0.5 * (sx[(int) d[1] - 1] + sx[(int) (d[1] + 0.5) - 1]);
  stats[2] = 0.5 * (sx[(int) d[2] - 1] + sx[(int) (d[2] + 0.5) - 1]);
  stats[3] = 0.5 * (sx[(int) d[3] - 1] + sx[(int) (d[3] + 0.5) - 1]);
  stats[4] = 0.5 * (sx[(int) d[4] - 1] + sx[(int) (d[4] + 0.5) - 1]);

  float iqr = stats[3] - stats[1];

  // Iterate through the values in sx and look for those that fall outside
  // of the range.
  float delta = coef * iqr;
  for (i = 0; i < n; i++) {
    if (sx[i] < stats[1] - delta || sx[i] > stats[3] + delta) {
      outliers->outliers[outliers->n] = sx[i];
      outliers->n++;
    }
  }

  // Free up memory
  free(sx);

  // Return the list of outliers.
  return outliers;
}
