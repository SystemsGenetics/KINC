#ifndef _PEARSON_
#define _PEARSON_

#include "../similarity.h"

/**
 *
 */
class PearsonSimilarity: public PairWiseSimilarity {
  public:
    PearsonSimilarity(PairWiseSet * pws, int * samples, int min_obs);
    ~PearsonSimilarity();

    void run();
};

#endif
