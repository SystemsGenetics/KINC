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

// SSYEV prototype
extern "C" void ssyev_(char* jobz, char* uplo, int* n, float* a, int* lda,
                       float* w, float* work, int* lwork, int* info);


int do_threshold(int argc, char *argv[]);

void print_threshold_usage();


/**
 * A base class for thresholding methods to inherit from.
 *
 * This class constructor will read arguments from the command-line, and
 * provides the getters and setters for commonly needed class members.
 */
class ThresholdMethod {

  protected:
    // The expression matrix object.
    EMatrix * ematrix;
    // The directory where the expression matrix is found
    char * input_dir;
    // Specifies the correlation method that was used: pc, mi, sc
    char method[10];


    // DATA FILTERS FOR CLUSTERED SIMILARITY DATA
    // ------------------------------------------
    // The maximum number of missing values in the comparision.
    int max_missing;
    // The minimum number of samples in a cluster.
    int min_cluster_size;


  public:
    ThresholdMethod(int argc, char *argv[]);
    ~ThresholdMethod();

    // GETTERS
    // -------
    char * getCorMethod() { return method; }
    int getMaxMissing() { return max_missing; }
    int getMinClusterSize() { return min_cluster_size; }

    // TO BE IMPLEMENTED BY THE CHILD CLASS
    // ------------------------------------
    double findThreshold();
};

/**
 * Implements Random Matrix Theory (RMT) Thresholding
 */

class RMTThreshold : public ThresholdMethod {

  private:
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

    double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size);
    double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace);
    // Calculates the eigenvalues of the given matrix.
    float * calculateEigen(float * smatrix, int size);
    //
    double * unfolding(float * e, int size, int m);
    // Removes duplicate eigenvalues from an array of eigenvalues.
    float * degenerate(float* eigens, int size, int* newSize);

    float * read_similarity_matrix_bin_file(float th, int * size);
    float * read_similarity_matrix_cluster_file(float th, int * size);

  public:
    RMTThreshold(int argc, char *argv[]);
    ~RMTThreshold();

    double findThreshold();

};

#endif

