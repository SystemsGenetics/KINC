#include "PairWiseSimilarity.h"

/**
 * Constructor
 *
 * @param const char * method
 *   The type of similarity test to perform. Valid values include 'pc' for
 *   Pearson's correlation, 'sc' for Spearman's correlation and 'mi' for
 *   Mutual Information.
 * @param PairWiseSet *pws
 *   The PairWiseSet which provides information about the two genes being
 *   compared.
 * @param int min_obs
 *   Optional.  The number of observations required to perform the
 *   similarity test. If the number of samples is less than this value the
 *   test will not be performed.  The default is 30 samples.
 */
PairWiseSimilarity::PairWiseSimilarity(PairWiseSet *pws, int min_obs = 30) {
  this->pws = pws;
  this->score = NAN;
  this->samples = NULL;
  this->min_obs = min_obs;

  init();
}
/**
 * Constructor
 *
 * @param const char * method
 *   The type of similarity test to perform. Valid values include 'pc' for
 *   Pearson's correlation, 'sc' for Spearman's correlation and 'mi' for
 *   Mutual Information.
 * @param PairWiseSet *pws
 *   The PairWiseSet which provides information about the two genes being
 *   compared.
 * @param int min_obs
 *   Optional.  The number of observations required to perform the
 *   similarity test. If the number of samples is less than this value the
 *   test will not be performed.  The default is 30 samples.
 * @param int * samples
 *   Optional. The set of samples to use when performing the pair-wise
 *   comparision.  This is an array of zero's and ones that indicates if the
 *   samples should be included in the test.  It must be the same size as the
 *   n_orig value of the pws argument. If this argument is not provided then
 *   the samples array that is part of the pws argument is used.
 */
PairWiseSimilarity::PairWiseSimilarity(PairWiseSet *pws, int min_obs, int * samples) {
  this->pws = pws;
  this->score = NAN;
  this->samples = samples;
  this->min_obs = min_obs;

  init();
}

/**
 *
 */
void PairWiseSimilarity::init() {
  this->a = (double *) malloc(sizeof(double) * pws->n_orig);
  this->b = (double *) malloc(sizeof(double) * pws->n_orig);
  this->n = 0;

  // If a samples list is provided then use that list.
  if (samples) {
    for (int i = 0; i < pws->n_orig; i++) {
      // If the sample is a 1 then include it in the a & b arrays that will be
      // used for testing. Also, make sure the pws sample is not missing.
      // We don't want the caller to accidentally include a missing value.
      if (samples[i] == 1 && pws->samples[i] == 1) {
        a[this->n] = pws->x_orig[i];
        b[this->n] = pws->y_orig[i];
        this->n++;
      }
    }
  }
  // If a samples list is not provided then use the clean list that
  // comes with the pws.
  else {
    for (int i = 0; i < pws->n_clean; i++) {
      a[this->n] = pws->x_clean[i];
      b[this->n] = pws->y_clean[i];
    }
    this->n = pws->n_clean;
  }
}

/**
 * Desctructor
 */
PairWiseSimilarity::~PairWiseSimilarity() {
  free(a);
  free(b);
}
