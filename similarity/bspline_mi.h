#ifndef _BSPLINE_MI_
#define _BSPLINE_MI_

#include <math.h>
#include "../similarity.h"

/**
 *
 */
class MISimilarity: public PairWiseSimilarity {
  private:
    double mi_bins;
    double mi_degree;
  public:
    MISimilarity(PairWiseSet * pws, int * samples, int min_obs, double mi_bins, double mi_degree);
    ~MISimilarity();

    void run();
};

#endif
