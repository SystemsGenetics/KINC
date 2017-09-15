#ifndef _MISIMILARITY_
#define _MISIMILARITY_

#include <math.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_bspline.h>
#include "PairWiseSimilarity.h"

/**
 *
 */
class MISimilarity: public PairWiseSimilarity {

  private:
    float mi_bins;
    float mi_degree;

    float calculateBSplineMI(float *v1, float *v2, int n, int m, int k,
        float xmin, float ymin, float xmax, float ymax);

  public:
    MISimilarity(PairWiseSet * pws, int min_obs, int * samples);
    MISimilarity(PairWiseSet * pws, int min_obs, float mi_bins, float mi_degree);
    MISimilarity(PairWiseSet * pws, int min_obs, int * samples, float mi_bins, float mi_degree);
    ~MISimilarity();

    void run();
};

#endif
