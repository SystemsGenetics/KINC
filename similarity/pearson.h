#ifndef _PEARSON_
#define _PEARSON_

#include "../similarity.h"

/**
 *
 */
class PearsonSimilarty: public PairWiseSimilarity {
  public:
    PearsonSimilarty(PairWiseSet * pws, int min_obs);
    ~PearsonSimilarty();

    void run();
};

#endif
