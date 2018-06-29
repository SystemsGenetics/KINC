#ifndef PAIRWISE_LINALG_H
#define PAIRWISE_LINALG_H
#include <ace/core/core.h>

namespace Pairwise
{
   typedef union {
      float s[2];
   } Vector2;

   void vectorInitZero(Vector2& a);
   void vectorAdd(Vector2& a, const Vector2& b);
   void vectorAdd(Vector2& a, float c, const Vector2& b);
   void vectorSubtract(Vector2& a, const Vector2& b);
   void vectorScale(Vector2& a, float c);
   float vectorDot(const Vector2& a, const Vector2& b);
   float vectorDiffNorm(const Vector2& a, const Vector2& b);

   typedef union {
      float s[4];
   } Matrix2x2;

   void matrixInitIdentity(Matrix2x2& M);
   void matrixInitZero(Matrix2x2& M);
   void matrixAdd(Matrix2x2& A, float c, const Matrix2x2& B);
   void matrixScale(Matrix2x2& A, float c);
   void matrixInverse(const Matrix2x2& A, Matrix2x2& B, float *p_det);
   void matrixProduct(const Matrix2x2& A, const Vector2& x, Vector2& b);
   void matrixOuterProduct(const Vector2& a, const Vector2& b, Matrix2x2& C);
}

#endif
