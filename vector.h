#ifndef _VECTOR_
#define _VECTOR_

#include <stdlib.h>
#include <math.h>

void quickSortF(float* l, int size);
void quickSortD(double* l, int size);

void swapD(double* l, int idx1, int idx2);
void swapF(float* l, int idx1, int idx2);


void remove_missing_paired(double *a, double *b, int n, double *a2, double *b2, int *n2);
#endif
