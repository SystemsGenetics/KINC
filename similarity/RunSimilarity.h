#ifndef _SIMILARITY_
#define _SIMILARITY_

#include <getopt.h>
#include <sys/stat.h>
#include "SpearmanSimilarity.h"
#include "PearsonSimilarity.h"
#include "MISimilarity.h"
#include "./clustering/MixtureModelClustering.h"
#include "../general/misc.h"

// a global variable for the number of rows in each output file
#define ROWS_PER_OUTPUT_FILE 10000
// the number of bins in the correlation value histogram
#define HIST_BINS 100

class RunSimilarity {

  private:
    // The expression matrix object.
    EMatrix * ematrix;
    // Specifies the method: sc, pc, mi.
    char * method;
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // An array holding the histogram of similarity scores.
    int * histogram;

    // Variables for the expression matrix
    // -----------------------------------
    // Indicates if the expression matrix has headers.
    int headers;
    // The input file name
    char *infilename;
    // The number of rows in the input ematrix file (including the header)
    int rows;
    // The number of cols in the input ematrix file.
    int cols;
    // Indicates if missing values should be ignored in the EMatrix file.
    int omit_na;
    // Specifies the value that represents a missing value.
    char *na_val;
    // Specifies the transformation function: log2, none.
    char func[10];

    // Variables for clustering
    // ------------------------
    // Indicates the clustering method to use.
    char * clustering;
    // The total number of jobs that will be run at once.
    int num_jobs;
    // The index of this job within the total jobs.  Must be
    // between 1 and num_jobs (no zero index).
    int job_index;

    // Variables for mutual information
    // --------------------------------
    // The number of bins for the B-spline estimate of MI.
    int mi_bins;
    // The degree of the B-spline function.
    int mi_degree;

    // Variables for Mixture Models
    // -----------------------------
    // The criterion model. E.g.  BIC, ICL, NEC, CV, DCV.
    char criterion[4];
    // The maximum number of clusters to allow per comparision.
    int max_clusters;

    // Variables for mean shift clustering.
    // --------------------------------
    double msc_bw1;
    double msc_bw2;

    void writeHistogram();
    // Calcualtes pair-wise similarity score the traditional way.
    void executeTraditional();

  public:
    RunSimilarity(int argc, char *argv[]);
    ~RunSimilarity();
    void execute();
    static void printUsage();

};

#endif
