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
#include <math.h>
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
    float thresholdStart;
    float thresholdStep;
    float chiSoughtValue;
    int minEigenVectorSize;
    float finalTH;
    float finalChi;
    float minTH;
    float minChi;
    float maxChi;

    float nnsdHistogramBin;
    float chiSquareTestThreshold;
    int minUnfoldingPace;
    int maxUnfoldingPace;

    float getNNSDChiSquare(float* eigens, int size);
    float getNNSDPaceChiSquare(float* eigens, int size, float bin, int pace);
    // Calculates the eigenvalues of the given matrix.
    float * calculateEigen(float * smatrix, int size);
    //
    float * unfolding(float * e, int size, int m);
    float * unfoldingByCDF(float * e, int size, int m);
    // Removes duplicate eigenvalues from an array of eigenvalues.
    float * degenerate(float* eigens, int size, int* newSize);

    float * read_similarity_matrix_bin_file(float th, int * size);
    float * read_similarity_matrix_cluster_file(float th, int * size);

  public:
    RMTThreshold(EMatrix * ematrix, char ** method, int num_methods, char * th_method,
        float thresholdStart, float thresholdStep, float chiSoughtValue,
        char * clustering, int min_cluster_size, int max_missing, int max_modes,
        float min_range);
    ~RMTThreshold();

    float findThreshold();

};

#endif
