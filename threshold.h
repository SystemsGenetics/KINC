#ifndef _THRESHOLD_
#define _THRESHOLD_

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <dirent.h>
#include "vector.h"
#include "similarity.h"


typedef struct{

  int perf;           // indicates if performance monitoring should be enabled
  int rows;           // the number of rows in the expression matrix
  int cols;           // the number of columns in the expression matrix
  char const * infilename;   // the input file name
  char method[10];    // specifies the method: cor, mi
  int headers;

  // Set some global params for the RMM stepping.
  float thresholdStart; // the threshold to start generating cutMatrices with (actually, +thresholdStep)
  float thresholdStep;  // the threshold step for each iteration of generation
  float chiSoughtValue; // the chiValue the while loop will end on

  double nnsdHistogramBin;
  double chiSquareTestThreshold;
  int minUnfoldingPace;
  int maxUnfoldingPace;
  int mimiModuleSize;
  double edHistogramBin; // Eigenvalue Histogram Bin size

  char fileprefix[1024]; // the input filename without the prefix
  char const * inputDir;

  int numGenes;          // the number of genes, n, in the nxn similarity matrix
  int numLinesPerFile;   // the number of lines per bin file

  // Holds an array flagging which genes meet the given threshold in the matrix.
  int * UsedFlag;

  // An array for indicating the index of the gene in the cut matrix.
  int * cutM_index;

  // The minimum size of the cut matrix.  If the cut matrix is smaller than
  // this size, a test will no be performed.
  int min_size;

} RMTParameters;


// SSYEV prototype
extern void ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda,
                    float* w, float* work, int* lwork, int* info );

/**
 * Function Prototypes
 */
int do_threshold(int argc, char *argv[]);

int find_threshold(RMTParameters *params);

float * read_similarity_matrix_bin_file(float th, int * size, RMTParameters *params);
float * read_similarity_matrix_cluster_file(float th, int * size, RMTParameters *params);

float* calculateEigen(float* mat, int size);

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, RMTParameters *params);

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace, RMTParameters *params);

void print_threshold_usage();
#endif

