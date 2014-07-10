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


typedef struct{

  int perf;           // indicates if performance monitoring should be enabled
  char *infilename;   // the input file name
  char method[10];    // specifies the method: cor, mi

  //set some global params for the RMM stepping
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
  char* inputDir;

  int numGenes;          // the number of genes, n, in the nxn similarity matrix
  int numLinesPerFile;   // the number of lines per bin file

  int* UsedFlag;         // holds an array flagging use of indicies in the matrix

  int* index1;           // had to rename because index was already taken

} RMTParameters;


// SSYEV prototype
extern void ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda,
                    float* w, float* work, int* lwork, int* info );

/**
 * Function Prototypes
 */
int do_threshold(int argc, char *argv[]);

int find_threshold(RMTParameters params);

float * read_similarity_matrix(float th, int * size, RMTParameters params);

void quickSortF(float* l, int size);

float* calculateEigen(float* mat, int size);

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, double bin, int minPace, int maxPace);

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace);

void print_threshold_usage();
#endif

