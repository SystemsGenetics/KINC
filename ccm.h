#ifndef _CCM_
#define _CCM_

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


typedef struct{

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

} CCMParameters;

/**
 * Constants
 */

// a global variable for the number of rows in each output file
static int ROWS_PER_OUTPUT_FILE = 10000;

// the number of bins in the correlation value histogram
static int HIST_BINS = 100;


// function prototypes
double ** load_ematrix(CCMParameters params);
int * init_histogram(CCMParameters params);
void print_histogram(CCMParameters params, int * histogram);
void calculate_MI(CCMParameters params, double ** data, int * histogram);
double calculateBSplineMI(double *v1, double *v2, int n, int m, int k, double xmin, double ymin, double xmax, double ymax);
void calculate_pearson(CCMParameters params, double ** data, int * histogram);
void print_usage();

#endif
