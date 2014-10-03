#include "vector.h"

/*
 * @param double* l
 * @param int idx1
 * @param int idx2
 */

void swapD(double* l, int idx1, int idx2) {
  double temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

/*
 * @param float* l
 * @param int idx1
 * @param int idx2
 */

void swapF(float* l, int idx1, int idx2) {
  float temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

/*
 * @param double* l
 * @param int size
 */

void quickSortD(double* l, int size){
  if (size <= 1) {
    return;
  }
  int pivIdx = (int) size / 1.618; // golden ratio
  double pivot = l[pivIdx];
  swapD(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for (i = 0; i < size - 1; i++) {
    if(l[i] < pivot){
      swapD(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapD(l, size-1, leftPlace);
  quickSortD(l, leftPlace);
  quickSortD(&l[leftPlace + 1], size - leftPlace - 1);
  return;
}

/*
 * @param float* l
 * @param int size
 */

void quickSortF(float* l, int size){
  if (size <= 1) {
    return;
  }
  int pivIdx = (int) size / 1.618; //golden ratio
  float pivot = l[pivIdx];
  swapF(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for (i = 0; i < size - 1; i++) {
    if(l[i] < pivot){
      swapF(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapF(l, size - 1, leftPlace);
  quickSortF(l,leftPlace);
  quickSortF(&l[leftPlace + 1], size - leftPlace - 1);
  return;
}
