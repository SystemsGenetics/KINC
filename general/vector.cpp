#include "vector.h"


/**
 * Returns an array indicating the order that the elements should be sorted.
 *
 * @param float *z
 * @param int n
 */
int * orderArray(float *z, int n) {
  int i, j, num_less;

  int * order = (int *) malloc(sizeof(int) * n);
  int used[n];

  for (i = 0; i < n; i++) {
    used[i] = 0;
  }
  for (i = 0; i < n; i++) {
    num_less = 0;
    for (j = 0; j < n; j++) {
      if (i == j) {
        continue;
      }
      if (z[i] > z[j]) {
        num_less++;
      }
    }
    while (used[num_less] == 1) {
      num_less++;
    }
    order[i] = num_less;
    used[num_less] = 1;
  }
  return order;
}
/*
 * @param float* l
 * @param int idx1
 * @param int idx2
 */

void swapI(int* l, int idx1, int idx2) {
  int temp = l[idx1];
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
 * @param int* l
 * @param int size
 */

void quickSortI(int* l, int size){
  if (size <= 1) {
    return;
  }
  int pivIdx = (int) size / 1.618; // golden ratio
  int pivot = l[pivIdx];
  swapI(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for (i = 0; i < size - 1; i++) {
    if(l[i] < pivot){
      swapI(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapI(l, size-1, leftPlace);
  quickSortI(l, leftPlace);
  quickSortI(&l[leftPlace + 1], size - leftPlace - 1);
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

/**
 * Removes values from two arrays when at least one is missing.
 *
 * @param float *a
 *   Vector 1
 * @param float *b
 *   Vector 2
 * @param int n
 *   The size of vectors a and b.
 * @param float *a2
 *   Vector 1 but with missing values removed. Must be initialized to
 *   at least the size of n.
 * @param float *b2
 *   Vector 2 but with missing values removed. Must be initialized to
 *   at least the size of n.
 * @param n2
 *   The new size of the vectors a2 band b2
 * @param int * kept
 *   An array of ones and zeros indicating if the element in a or b were
 *   included in a2 and b2.
 */
void remove_missing_paired(float *a, float *b, int n, float *a2, float *b2, int *n2, int *kept) {
  int i;
  *n2 = 0;
  for (i = 0; i < n; i++) {
    // if either of these elements is missing then don't include the
    // elements from this sample
    if (isnan(a[i]) || isnan(b[i]) || isinf(a[i]) || isinf(b[i])) {
      kept[i] = 0;
      continue;
    }
    a2[*n2] = a[i];
    b2[*n2] = b[i];
    *n2 = *n2 + 1;
    kept[i] = 1;
  }
}
