#ifndef _PAIRWISESET_
#define _PAIRWISESET_

#include "../ematrix/EMatrix.h"
#include "../stats/outlier.h"

/**
 * A class that holds expression data for two genes/probesets.
 *
 */
class PairWiseSet {

  private:
    void clean();

  public:
    // The indexes into the EMatrix for the two genes being compared.
    int gene1;
    int gene2;
    // The original x and y data arrays and their size.
    double *x_orig;
    double *y_orig;
    int n_orig;
    // The x and y data arrays after NAs have been removed and their size.
    double *x_clean;
    double *y_clean;
    int n_clean;
    // The samples array of zeros and ones. Where zero indicates the sample
    // was removed and not included in the comparision and one indicates it
    // was preserved and used in the comparision.
    int * samples;
    // The threshold for expression values.
    double threshold;

  public:
    PairWiseSet(EMatrix * ematrix, int i, int j);
    PairWiseSet(EMatrix * ematrix, int i, int j, double th);
    PairWiseSet(double *a, double *b, int n, int i, int j);
    PairWiseSet(double *a, double *b, int n, int i, int j, double th);
    ~PairWiseSet();

    // Detects and removes outliers from the input vectors.
    void maskOutliers();

};

#endif
