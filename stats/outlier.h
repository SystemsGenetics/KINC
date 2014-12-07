#ifndef _OUTLIER_
#define _OUTLIER_

#include "../error.h"
#include "../vector.h"


typedef struct {
  // The list of values that are considered outliers.
  double * outliers;
  // The number of outliers
  int n;
} Outliers;


Outliers outliers_iqr(double * x, int n, double coef);

#endif
