#ifndef _SPEARMAN_
#define _SPEARMAN_

#include "../similarity.h"


/**
 * Class for Spearman Correlation similarity.
 */
class SpearmanSimilarty: public PairWiseSimilarity {
  public:
    SpearmanSimilarty(PairWiseSet * pws, int min_obs);
    ~SpearmanSimilarty();

    void run();
};

/**
 *
 */
void calculate_spearman(CCMParameters params, double ** data, int * histogram);

#endif
