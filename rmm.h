#ifndef _RANDOMMATRIX_
#define _RANDOMMATRIX_

typedef struct{

  double nnsdHistogramBin;
  double chiSquareTestThreshold;
  int minUnfoldingPace;
  int maxUnfoldingPace;
  int mimiModuleSize;
  double edHistogramBin; // Eigenvalue Histogram Bin size

} RMTParameters;

int determinePearsonCorrelationThreshold_LargeMatrix();

float * readPearsonCorrelationMatrix(float th, int * size);

FILE * fileOpenHelper(char* extension);

// prototype declarations
void quickSortF(float* l, int size);

float* calculateEigen(float* mat, int size);

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, double bin, int minPace, int maxPace);

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace);

#endif

