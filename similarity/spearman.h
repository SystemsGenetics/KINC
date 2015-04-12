#ifndef _SPEARMAN_
#define _SPEARMAN_

#include "../similarity.h"


/**
 * Class for Spearman Correlation similarity.
 */
class SpearmanSimilarity: public PairWiseSimilarity {
  public:
    SpearmanSimilarity(PairWiseSet * pws, int * samples, int min_obs);
    ~SpearmanSimilarity();

    void run();
};

/**
 *
 */
void calculate_spearman(CCMParameters params, double ** data, int * histogram);

#endif
