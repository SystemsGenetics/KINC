#ifndef _BSPLINE_MI_
#define _BSPLINE_MI_

#include "../similarity.h"

void calculate_MI(CCMParameters params, double ** data, int * histogram);
double calculateBSplineMI(double *v1, double *v2, int n, int m, int k, double xmin, double ymin, double xmax, double ymax);

#endif
