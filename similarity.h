#ifndef _SIMILARITY_
#define _SIMILARITY_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_bspline.h>

#include "ematrix.h"

/**
 * A class that holds expression data for two genes/probesets.
 *
 */
class PairWiseSet {
  // The PairWiseSimilarity class should have access to the protected
  // members of this class.
  friend class PairWiseSimilarity;
  friend class MISimilarity;
  friend class PearsonSimilarity;
  friend class SpearmanSimilarity;
  friend class PairWiseClusterWriter;

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

  public:
    PairWiseSet(EMatrix * ematrix, int i, int j);
    PairWiseSet(double *a, double *b, int n, int i, int j);
    ~PairWiseSet();
};

/**
 * A base class for similiarity functions.
 *
 * The pair-wise comparision is performed on a single PairWiseSet.
 */
class PairWiseSimilarity {
  friend class PairWiseClusterWriter;
  protected:
    // The PairWiseSet containing the two samples for comparision.
    PairWiseSet *pws;
    // The method name (e.g. pc, mi, sc).
    char * method;
    // The final similarity score.
    double score;
    // An array of 1's and zeros indicating which samples should be
    // included in the pair-wise comparision.
    int * samples;
    // The minimum number of observations required.
    int min_obs;

    // The expression arrays with non included samples removed and the
    // size of these arrays.
    double *a, *b;
    int n;

  public:
    PairWiseSimilarity(const char * method, PairWiseSet *pws, int * samples, int min_obs);
    ~PairWiseSimilarity();

    // Executes the pair-wise similiarity function. This should
    // be implemented by the descendent class.
    void run() {}
};


typedef struct {

  int perf;          // indicates if performance monitoring should be enabled
  int omit_na;       // indicates if missing values should be ignored
  int headers;

  // Command-line option values
  int rows;           // the number of rows in the expression matrix
  int cols;           // the number of columns in the expression matrix
  char *infilename;   // the input file name
  char *na_val;       // specifies the value that represents a missing value
  char func[10];      // specifies the transformation function: log2, none
  char method[10];    // specifies the method: cor, mi
  int min_obs;        // the minimum number of observations to calculate correlation


  // other global variables
  int do_log10;          // set to 1 to perform log10 transformation
  int do_log2;           // set to 1 to perform log2 transformation
  int do_log;            // set to 1 to perform log transformation
  char fileprefix[1024]; // the input filename without the prefix

  // variables for mutual information
  int mi_bins;        // the number of bins for the B-spline estimate of MI
  int mi_degree;      // the degree of the B-spline function

  // variables for mean shift clustering
  double msc_bw1;
  double msc_bw2;

} CCMParameters;

// a global variable for the number of rows in each output file
#define ROWS_PER_OUTPUT_FILE 10000

// the number of bins in the correlation value histogram
#define HIST_BINS 100

// function prototypes

int do_similarity(int argc, char *argv[]);

void print_similarity_usage();

int * init_histogram(CCMParameters params);

void print_histogram(CCMParameters params, int * histogram);

void calculate_pearson(CCMParameters params, double ** data, int * histogram);

void calculate_spearman(CCMParameters params, double ** data, int * histogram);

void calculate_MI(CCMParameters params, double ** data, int * histogram);

double calculateBSplineMI(double *v1, double *v2, int n, int m, int k, double xmin, double ymin, double xmax, double ymax);


#endif
