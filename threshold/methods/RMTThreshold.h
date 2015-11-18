#ifndef _RMTTHRESHOLD_
#define _RMTTHRESHOLD_


// GNU Scientific Library headers.
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <sys/stat.h>
#include <dirent.h>
#include <regex.h>
#include "../../general/vector.h"


#include "ThresholdMethod.h"

// SSYEV prototype
extern "C" void ssyev_(char* jobz, char* uplo, int* n, float* a, int* lda,
                       float* w, float* work, int* lwork, int* info);

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
    RMTThreshold(EMatrix * ematrix, char ** method, int num_methods, char * th_method,
        double thresholdStart, double thresholdStep, double chiSoughtValue,
        char * clustering, int min_cluster_size, int max_missing, int max_modes);
    ~RMTThreshold();

    double findThreshold();

};

#endif
