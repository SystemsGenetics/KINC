#ifndef _SPEARMAN_
#define _SPEARMAN_

#include <gsl/gsl_statistics.h>
#include "PairWiseSimilarity.h"

/**
 * Class for Spearman Correlation similarity.
 */
class SpearmanSimilarity : public PairWiseSimilarity {
  public:
    SpearmanSimilarity(PairWiseSet * pws, int min_obs);
    SpearmanSimilarity(PairWiseSet * pws, int min_obs, int * samples);
    ~SpearmanSimilarity();

    void run();
};

#endif
