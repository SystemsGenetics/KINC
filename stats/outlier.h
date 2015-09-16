#ifndef _OUTLIER_
#define _OUTLIER_

#include "../general/error.h"
#include "../general/vector.h"


typedef struct {
  // The list of values that are considered outliers.
  double * outliers;
  // The number of outliers
  int n;
} Outliers;


Outliers * outliers_iqr(double * x, int n, double coef);

#endif
