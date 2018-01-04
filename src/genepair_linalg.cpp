#include "genepair_linalg.h"



namespace GenePair {






#define ELEM(M, n, i, j) ((M)[(i) * (n) + (j)])






void vectorInitZero(float *a)
{
   a[0] = 0;
   a[1] = 0;
}






void vectorCopy(float *a, const float *b)
{
   a[0] = b[0];
   a[1] = b[1];
}






void vectorAdd(float *a, const float *b)
{
   a[0] += b[0];
   a[1] += b[1];
}






void vectorAdd(float *a, float c, const float *b)
{
   a[0] += c * b[0];
   a[1] += c * b[1];
}






void vectorSubtract(float *a, const float *b)
{
   a[0] -= b[0];
   a[1] -= b[1];
}






void vectorScale(float *a, float c)
{
   a[0] *= c;
   a[1] *= c;
}






float vectorDot(const float *a, const float *b)
{
   return a[0] * b[0] + a[1] * b[1];
}






float vectorDiffNorm(const float *a, const float *b)
{
   float dist = 0;
   dist += (a[0] - b[0]) * (a[0] - b[0]);
   dist += (a[1] - b[1]) * (a[1] - b[1]);

   return sqrt(dist);
}






void matrixInitIdentity(float *M)
{
   const int N = 2;

   ELEM(M, N, 0, 0) = 1;
   ELEM(M, N, 0, 1) = 0;
   ELEM(M, N, 1, 0) = 0;
   ELEM(M, N, 1, 1) = 1;
}






void matrixInitZero(float *M)
{
   const int N = 2;

   ELEM(M, N, 0, 0) = 0;
   ELEM(M, N, 0, 1) = 0;
   ELEM(M, N, 1, 0) = 0;
   ELEM(M, N, 1, 1) = 0;
}






void matrixAdd(float *A, float c, const float *B)
{
   const int N = 2;

   ELEM(A, N, 0, 0) += c * ELEM(B, N, 0, 0);
   ELEM(A, N, 0, 1) += c * ELEM(B, N, 0, 1);
   ELEM(A, N, 1, 0) += c * ELEM(B, N, 1, 0);
   ELEM(A, N, 1, 1) += c * ELEM(B, N, 1, 1);
}






void matrixScale(float *A, float c)
{
   const int N = 2;

   ELEM(A, N, 0, 0) *= c;
   ELEM(A, N, 0, 1) *= c;
   ELEM(A, N, 1, 0) *= c;
   ELEM(A, N, 1, 1) *= c;
}






void matrixInverse(const float *A, float *B, float *p_det)
{
   const int N = 2;
   const float EPSILON = 1e-5;

   float det = ELEM(A, N, 0, 0) * ELEM(A, N, 1, 1) - ELEM(A, N, 0, 1) * ELEM(A, N, 1, 0);

   if ( fabs(det) <= EPSILON )
   {
      throw std::runtime_error("singular matrix");
   }

   ELEM(B, N, 0, 0) = +ELEM(A, N, 1, 1) / det;
   ELEM(B, N, 0, 1) = -ELEM(A, N, 0, 1) / det;
   ELEM(B, N, 1, 0) = -ELEM(A, N, 1, 0) / det;
   ELEM(B, N, 1, 1) = +ELEM(A, N, 0, 0) / det;

   *p_det = det;
}






void matrixProduct(const float *A, const float *x, float *b)
{
   const int N = 2;

   b[0] = ELEM(A, N, 0, 0) * x[0] + ELEM(A, N, 0, 1) * x[1];
   b[1] = ELEM(A, N, 1, 0) * x[0] + ELEM(A, N, 1, 1) * x[1];
}






void matrixOuterProduct(const float *a, const float *b, float *C)
{
   const int N = 2;

   ELEM(C, N, 0, 0) = a[0] * b[0];
   ELEM(C, N, 0, 1) = a[0] * b[1];
   ELEM(C, N, 1, 0) = a[1] * b[0];
   ELEM(C, N, 1, 1) = a[1] * b[1];
}






}
