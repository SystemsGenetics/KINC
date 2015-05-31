#ifndef _SIMILARITY_
#define _SIMILARITY_

#include <getopt.h>
#include <sys/stat.h>
#include "SpearmanSimilarity.h"
#include "PearsonSimilarity.h"
#include "MISimilarity.h"

// a global variable for the number of rows in each output file
#define ROWS_PER_OUTPUT_FILE 10000
// the number of bins in the correlation value histogram
#define HIST_BINS 100

void print_similarity_usage();

class Similarity {

  private:
    // The expression matrix object.
    EMatrix * ematrix;
    // Specifies the method: sc, pc, mi.
    char method[10];
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // An array holding the histogram of similarity scores.
    int * histogram;


    // Variables for mutual information
    // --------------------------------
    // The number of bins for the B-spline estimate of MI.
    int mi_bins;
    // The degree of the B-spline function.
    int mi_degree;

    // Variables for mean shift clustering.
    // --------------------------------
//    double msc_bw1;
//    double msc_bw2;

    void writeHistogram();

  public:
    Similarity(int argc, char *argv[]);
    ~Similarity();
    void run();

};

#endif
