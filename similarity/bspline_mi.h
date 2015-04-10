#ifndef _BSPLINE_MI_
#define _BSPLINE_MI_

#include <math.h>
#include "../similarity.h"

/**
 *
 */
class MISimilarty: public PairWiseSimilarity {
  private:
    double mi_bins;
    double mi_degree;
  public:
    MISimilarty(PairWiseSet * pws, int min_obs, double mi_bins, double mi_degree);
    ~MISimilarty();

    void run();
};

#endif
