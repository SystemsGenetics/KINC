#ifndef _PAIRWISESIMILARITY_
#define _PAIRWISESIMILARITY_

#include <getopt.h>
#include "../PairWiseSet.h"

/**
 * A base class for similiarity functions.
 *
 * The pair-wise comparision is performed on a single PairWiseSet.
 */
class PairWiseSimilarity {
  protected:
    // The PairWiseSet containing the two samples for comparision.
    PairWiseSet *pws;
    // The final similarity score.
    double score;
    // An array of 1's and zeros indicating which samples should be
    // included in the pair-wise comparision.
    int * samples;
    // The expression arrays with non included samples removed and the
    // size of these arrays.
    double *a, *b;
    // The number of samples in this pair-wise test.
    int n;
    // The minimum number of observations required to perform the comparision.
    int min_obs;

    // The type of similarity that was performed.
    char * type;

    // Called by the constructors to initialize the object.
    void init();

  public:
    PairWiseSimilarity(PairWiseSet *pws, int min_obs);
    PairWiseSimilarity(PairWiseSet *pws, int min_obs, int * samples);
    ~PairWiseSimilarity();

    // Get's the score
    double getScore() { return score; }

    // Executes the pair-wise similiarity function. This should
    // be implemented by the descendent class.
    void run() {}

};

#endif
