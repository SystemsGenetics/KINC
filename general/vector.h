#ifndef _VECTOR_
#define _VECTOR_

#include <stdlib.h>
#include <math.h>

int * orderArray(float *z, int n);

void quickSortF(float* l, int size);
void quickSortI(int* l, int size);

void swapF(float* l, int idx1, int idx2);
void swapI(int* l, int idx1, int idx2);

void remove_missing_paired(float *a, float *b, int n, float *a2, float *b2, int *n2, int *kept);

#endif
