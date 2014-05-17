/* v2: Replaces the GSL eigenvalue solver with the MAGMA solver */
/* v2-1: Branching from v2, attempting to use the ACML lapack
 *       eigenvector solver instead of MAGMA (or perhaps as an 
 *       intermediate step. */
/* 
 * v2-2: Uses MKL rather than ACML for eigen solve */
/* 
 * v2-2.1: Removed the sort call and eigen print inside of chiSquareTest */
/*
 * v2-3:   Uses a symmettric matrix eigen solver.  Use with RMM.v1-3 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

/* SSYEV prototype */
extern void ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda,
                    float* w, float* work, int* lwork, int* info );

void swapD(double* l, int idx1, int idx2){
  double temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

void swapF(float* l, int idx1, int idx2){
  float temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

void quickSortD(double* l, int size){
  if(size<=1) return;
  int pivIdx = (int) size/1.618;//golden ratio
  double pivot = l[pivIdx];
  swapD(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for(i=0;i<size-1;i++){
    if(l[i]<pivot){
      swapD(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapD(l, size-1, leftPlace);
  quickSortD(l,leftPlace);
  quickSortD(&l[leftPlace+1], size-leftPlace-1);
  return;
}

void quickSortF(float* l, int size){
  if(size<=1) return;
  int pivIdx = (int) size/1.618;//golden ratio
  float pivot = l[pivIdx];
  swapF(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for(i=0;i<size-1;i++){
    if(l[i]<pivot){
      swapF(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapF(l, size-1, leftPlace);
  quickSortF(l,leftPlace);
  quickSortF(&l[leftPlace+1], size-leftPlace-1);
  return;
}

/**
 * Calculates the eigenvalues of the given matrix.  This function is a wrapper
 * for the ssyev_ function of the LAPACK package.
 *
 * @param float *mat
 *   A pointer to an array of floating point numbers representing the
 *   square n x n correlation matrix.
 * @param int size
 *   The size, n, of the cut n x n matrix.
 *
 * @return
 *   A pointer to an array of floating point numbers representing the 
 *   array of eigenvalues
 */
float* calculateEigen(float* mat, int size){

  char jobz = 'N';      // N means to not don't compute eigenvectors, just eigenvalues
  char uplo = 'U';      // U means the upper matrix is stored
  float * W;            // the array where eignvalues are stored
  float * work;         // a working array. This will be 5 times the size of the final array
  int lwork = 5 * size; // the size of the work array
  int rc;               // indicates the success of the ssyev_ function

  // allocate the arrays 
  W    = (float *) malloc(sizeof(float) * size);
  work = (float *) malloc(sizeof(float) * 5 * size);

  ssyev_(&jobz, &uplo , &size, mat, &size, W, work, &lwork, &rc);

  // report any errors
  if (rc < 0) {
    printf("\nERROR: During eigenvalue calculation, the %d argument had an illegal value. Continuing anyway...\n", rc);
  }
  else if (rc > 0) {
    printf("\nERROR: The eigenvalue algorithm failed to converge; %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero. Continuing anyway...\n", rc);
  }
  free(work);
  return W;
}

//returned array will always be sorted and of length size-1
double* unfolding(float* e, int size, int m){
  int count=1, i,j=0;//count equals 1 initially because of 2 lines following loop which propogates the arrays
  for(i=0; i<size-m; i+=m) count++;
  double* oX = (double*) malloc(sizeof(double)*count);
  double* oY = (double*) malloc(sizeof(double)*count);
  for(i=0; i<size-m; i+=m){
    oX[j]=e[i];
    oY[j]= (i+1.0)/(double)size;
    j++;
  }
  oX[count-1] = e[size-1];
  oY[count-1] = 1;

  for(i=1;i<count;i++){
    if(!(oX[i-1]<oX[i])){
      printf("\nat postion %d a problem exists\n", i);
      printf("oX[i-1]=%f whilst oX[i]=%f\n",oX[i-1],oX[i]);
    }
  }
  double* yy = (double*) malloc(sizeof(double)*(size));

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, count);//see gsl docs, chapter 27: cspline is a natural spline
  gsl_spline_init(spline, oX, oY, count);

  for(i=0;i<(size-2);i++){
    yy[i+1] = gsl_spline_eval(spline, e[i+1], acc);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  yy[0] = 0.0;
  yy[size-1] = 1.0;
  for(i=0;i<size-1;i++){
    yy[i] = size*(yy[i+1]-yy[i]);
  }
  quickSortD(yy, size-1);
  free(oX);
  free(oY);
  return yy;
}

float* degenerate(float* eigens, int size, int* newSize){
  int i, j=0, count=1;//because one flag is set before the loop
  for(i=0;i<size;i++){
    if(fabs(eigens[i]) < 0.000001){
      eigens[i] = 0.0;
    }
  }
  int* flags = (int*) malloc(sizeof(int)*size);
  memset(flags, 0, size*sizeof(int));
  float temp = eigens[0];
  flags[0]=1;
  for(i=1;i<size; i++){
    if(fabs(eigens[i]-temp) > 0.000001){
      count++;
      flags[i] = 1;
      temp = eigens[i];
    }
  }
  float* remDups = (float*) malloc(sizeof(float)*count);//remDups means "removed duplicates"
  for(i=0;i<size;i++){
    if(flags[i]==1){
      remDups[j] = eigens[i];
      j++;
    }
  }
  free(flags);
  *newSize = count;

  return remDups;
}


double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace){
  int newSize;
  float* newE;
  double* edif;
  newE = degenerate(eigens, size, &newSize);
  size = newSize;

  edif = unfolding(newE, size, pace);
  free(newE);
  size = size-1; //see note above unfolding function, will return an array of size-1
  int n = (int) (3.0/bin) + 1;
  double obj, expect, chi = 0;
  int i, j, count;
  for(i=0;i<n;i++){
    count = 0;
    for(j=0;j<size;j++){
      if(edif[j]>i*bin && edif[j] < (i+1)*bin) count++;
    }
    obj = (double) count;
    expect = (exp(-1*i*bin)-exp(-1*(i+1)*bin))*size;
    chi += (obj -expect)*(obj-expect)/expect;
  }
  free(edif);
  return chi;
}

//calls same name, 4 args instead of 5
double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, double bin, int minPace, int maxPace){
  double chiTest =0;
  int i = 0;
  int m;

  i=0;
  for(m = minPace; m<maxPace; m++){
    chiTest+=chiSquareTestUnfoldingNNSDWithPoisson4(eigens, size, bin, m);
    i++;
  }

  return chiTest/i;
}

