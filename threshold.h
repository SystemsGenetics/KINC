#ifndef _THRESHOLD_
#define _THRESHOLD_

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <dirent.h>

// GNU Scientific Library headers.
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

// KINC headers.
#include "vector.h"
#include "similarity.h"

/**
 * This class holds the arguments used as input to the threshold program.
 */
class RMTArgs {

  private:
    // Indicates if headers are present in the input EMatrix file.
    int headers;
    // The number of rows in the expression matrix.
    int rows;
    // The number of columns in the expression matrix.
    int cols;
    // The input file name
    char *infilename;
    // Specifies the correlation method: pc, mi, sc
    char method[10];
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // The input filename without the prefix.
    char * fileprefix;
    // The input directory where ematrix file is stored.
    char * inputDir;
    // The total number of jobs that will be run at once.
    int num_jobs;
    // The index of this job within the total jobs.  Must be
    // between 1 and num_jobs (no zero index)
    int job_index;

    // Filters
    int max_missing;
    int min_cluster_size;

    // Variables for RMT
    double thresholdStart;
    double thresholdStep;
    double chiSoughtValue;
    int minEigenVectorSize;
    double finalTH;
    double finalChi;
    double minTH;
    double minChi;
    double maxChi;

    double nnsdHistogramBin;
    double chiSquareTestThreshold;
    int minUnfoldingPace;
    int maxUnfoldingPace;
    int mimiModuleSize;
    double edHistogramBin; // Eigenvalue Histogram Bin size

  public:
    RMTArgs(int argc, char *argv[]);
    ~RMTArgs();

    // Getters
    int getHasHeaders() { return headers; }
    int getNumRows() { return rows; }
    int getNumCols() { return cols; }
    char * getInfileName() { return infilename; }
    char * getCorMethod() { return method; }
    char * getFilePrefix() { return fileprefix; }
    char * getInputDir() { return inputDir; }

    double getThresholdStart() { return thresholdStart; };
    double getThresholdStep() { return thresholdStep; };
    double getChiSoughtValue() { return chiSoughtValue; };
    int getMinEigenVectorSize() { return minEigenVectorSize; };
    double getFinalTH() { return finalTH; };
    double getFinalChi() { return finalChi; };
    double getMinTH() { return minTH; };
    double getMinChi() { return minChi; };
    double getMaxChi() { return maxChi; };
    double getNNSDHistogramBin() { return nnsdHistogramBin; }
    double getChiSquareTestThreshold() { return chiSquareTestThreshold; }
    int getMinUnfoldingPace() { return minUnfoldingPace; }
    int getMaxUnfoldingPace() { return maxUnfoldingPace; }
    int getMaxMissing() { return max_missing; }
    int getMinClusterSize() { return min_cluster_size; }

    int getMimiModuleSize() { return mimiModuleSize; }
    double getEdHistogramBin() { return edHistogramBin; }
};


// SSYEV prototype
extern "C" void ssyev_(char* jobz, char* uplo, int* n, float* a, int* lda,
                       float* w, float* work, int* lwork, int* info);

/**
 * Function Prototypes
 */
int do_threshold(int argc, char *argv[]);

int find_threshold(RMTArgs *params);

float * read_similarity_matrix_bin_file(float th, int * size, RMTArgs *params);
float * read_similarity_matrix_cluster_file(float th, int * size, RMTArgs *params);

float* calculateEigen(float* mat, int size);

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, RMTArgs *params);

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace, RMTArgs *params);

void print_threshold_usage();
#endif

