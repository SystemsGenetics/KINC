#ifndef _PEARSON_
#define _PEARSON_

#include <gsl/gsl_statistics.h>
#include "PairWiseSimilarity.h"

/**
 *
 */
class PearsonSimilarity: public PairWiseSimilarity {
  public:
    PearsonSimilarity(PairWiseSet * pws, int min_obs);
    PearsonSimilarity(PairWiseSet * pws, int min_obs, int * samples);
    ~PearsonSimilarity();

    void run();
};

#endif
