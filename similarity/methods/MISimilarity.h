#ifndef _MISIMILARITY_
#define _MISIMILARITY_

#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_bspline.h>
#include "PairWiseSimilarity.h"

/**
 *
 */
class MISimilarity: public PairWiseSimilarity {

  private:
    double mi_bins;
    double mi_degree;

    double calculateBSplineMI(double *v1, double *v2, int n, int m, int k,
        double xmin, double ymin, double xmax, double ymax);

  public:
    MISimilarity(PairWiseSet * pws, int min_obs, int * samples);
    MISimilarity(PairWiseSet * pws, int min_obs, double mi_bins, double mi_degree);
    MISimilarity(PairWiseSet * pws, int min_obs, int * samples, double mi_bins, double mi_degree);
    ~MISimilarity();

    void run();
};

#endif
