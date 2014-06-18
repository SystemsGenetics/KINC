#ifndef _RANDOMMATRIX_
#define _RANDOMMATRIX_

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
#include <setjmp.h>

typedef struct{

  double nnsdHistogramBin;
  double chiSquareTestThreshold;
  int minUnfoldingPace;
  int maxUnfoldingPace;
  int mimiModuleSize;
  double edHistogramBin; // Eigenvalue Histogram Bin size

} RMTParameters;

/**
 * Definitions for mimicing a try, catch block.
 */
#define TRY do{ jmp_buf ex_buf__; if( !setjmp(ex_buf__) ){
#define CATCH } else {
#define ETRY } }while(0)
#define THROW longjmp(ex_buf__, 1)


/**
 * Function Prototypes
 */
int determinePearsonCorrelationThreshold_LargeMatrix();

float * readPearsonCorrelationMatrix(float th, int * size);

FILE * fileOpenHelper(char* extension);

// prototype declarations
void quickSortF(float* l, int size);

float* calculateEigen(float* mat, int size);

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, double bin, int minPace, int maxPace);

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace);

#endif

