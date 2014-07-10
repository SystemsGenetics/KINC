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

/**
 * Description:
 * ------------
 * Calculates a Pearson correlation or Mutual Information matrix from an n x m
 * expression matrix where the columns are the samples and rows are
 * the genes or probesets. Cells represent expression levels. The first
 * column of the expression matrix should be the name of the gene or
 * probeset.
 *
 *
 * Change Log / History:
 * ----------
 * 6/2014 (by Stephen Ficklin)
 * 1) Added support for calculation of mutual information matrix
 * 2) Use getopt_long for --xxxx and -x style input arguments
 * 3) Re-organized the code
 *
 * 10/2014 (by Stephen Ficklin)
 * 1) Input files no longer require an implied .txt extesion but the full filename
 *    should be specified
 * 2) Removed the spaces in the output file name and renamed file.
 * 3) Incorporated GSL Pearson calculation function because the pre-calculating
 *    of sum of row and sum of squared of the row couldn't account for missing values
 *
 * 10/2011
 * Created by Scott Gibson at Clemson University under direction of
 * Dr. Melissa Smith in Collaboration with Alex Feltus, Feng Luo and Stephen
 * Ficklin
 */

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

// a global variable for the number of rows in each output file
static int ROWS_PER_OUTPUT_FILE = 10000;

// the number of bins in the correlation value histogram
static int HIST_BINS = 100;

// function prototypes

int do_similarity(int argc, char *argv[]);

double ** load_ematrix(CCMParameters params);

void print_similarity_usage();

int * init_histogram(CCMParameters params);

void print_histogram(CCMParameters params, int * histogram);

#endif
