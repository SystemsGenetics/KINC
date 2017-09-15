#ifndef _SIMILARITY_
#define _SIMILARITY_

#include <getopt.h>
#include <sys/stat.h>
#include <string.h>
#include "GeneFilter.h"
#include "./methods/SpearmanSimilarity.h"
#include "./methods/PearsonSimilarity.h"
#include "./methods/MISimilarity.h"
#include "./clustering/MixtureModelClustering.h"
#include "../general/misc.h"
#include "../stats/royston.h"

// a global variable for the number of rows in each output file
#define ROWS_PER_OUTPUT_FILE 10000
// the number of bins in the correlation value histogram
#define HIST_BINS 100

class RunSimilarity {

  private:
    // The expression matrix object.
    EMatrix * ematrix;
    // Specifies the methods: sc, pc, mi.
    char ** method;
    // Indicates the number of methods.
    int num_methods;
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // An array holding the histogram of similarity scores.
    int * histogram;
    // The threshold for expression values.
    float threshold;

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
    // An array that specifies the set of minimum thresholds to include.
    // All methods must have a minimum value, -1 should be used to exclude
    // a filter for that method.
    float * min_sim;

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
    float msc_bw1;
    float msc_bw2;

    // Variables to limit pair-wise comparisions to a specific set
    // of genes.
    // --------------------------------
    // An array of gene indexes that should be used for correlation. If this
    // list is set then only these genes will have a pair-wise similarity test
    // with every other gene.  If list2 is also set then pair-wise correlation
    // of genes in list1 will only occur with genes in list2.
    geneFilter set1;

    // A second array of gene indexes that should be used for correlation.
    // This list can only be set if the list1 is also set.
    geneFilter set2;

    void writeHistogram();
    // Calcualtes pair-wise similarity score the traditional way.
    void executeTraditional();
    void parseMethods(char * methods_str);
    void getGeneSet(geneFilter * set);
    void parseMinSim(char * minsim_str);

  public:
    RunSimilarity(int argc, char *argv[]);
    ~RunSimilarity();
    void execute();
    static void printUsage();

};

#endif
