#ifndef GENEPAIR_LINALG_H
#define GENEPAIR_LINALG_H
#include <ace/core/AceCore.h>

namespace GenePair
{
   void vectorInitZero(float *a);
   void vectorCopy(float *a, const float *b);
   void vectorAdd(float *a, const float *b);
   void vectorAdd(float *a, float c, const float *b);
   void vectorSubtract(float *a, const float *b);
   void vectorScale(float *a, float c);
   float vectorDot(const float *a, const float *b);
   float vectorDiffNorm(const float *a, const float *b);

   void matrixInitIdentity(float *M);
   void matrixInitZero(float *M);
   void matrixAdd(float *A, float c, const float *B);
   void matrixScale(float *A, float c);
   void matrixInverse(const float *A, float *B, float *p_det);
   void matrixProduct(const float *A, const float *x, float *b);
   void matrixOuterProduct(const float *a, const float *b, float *C);
}

#endif
