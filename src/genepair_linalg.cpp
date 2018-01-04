#include "genepair_linalg.h"



namespace GenePair {






inline const float& elem(const Matrix2x2& M, int n, int i, int j)
{
   return M.data[i * n + j];
}






inline float& elem(Matrix2x2& M, int n, int i, int j)
{
   return M.data[i * n + j];
}






void vectorInitZero(Vector2& a)
{
   a.data[0] = 0;
   a.data[1] = 0;
}






void vectorAdd(Vector2& a, const Vector2& b)
{
   a.data[0] += b.data[0];
   a.data[1] += b.data[1];
}






void vectorAdd(Vector2& a, float c, const Vector2& b)
{
   a.data[0] += c * b.data[0];
   a.data[1] += c * b.data[1];
}






void vectorSubtract(Vector2& a, const Vector2& b)
{
   a.data[0] -= b.data[0];
   a.data[1] -= b.data[1];
}






void vectorScale(Vector2& a, float c)
{
   a.data[0] *= c;
   a.data[1] *= c;
}






float vectorDot(const Vector2& a, const Vector2& b)
{
   return a.data[0] * b.data[0] + a.data[1] * b.data[1];
}






float vectorDiffNorm(const Vector2& a, const Vector2& b)
{
   float dist = 0;
   dist += (a.data[0] - b.data[0]) * (a.data[0] - b.data[0]);
   dist += (a.data[1] - b.data[1]) * (a.data[1] - b.data[1]);

   return sqrt(dist);
}






void matrixInitIdentity(Matrix2x2& M)
{
   const int N = 2;

   elem(M, N, 0, 0) = 1;
   elem(M, N, 0, 1) = 0;
   elem(M, N, 1, 0) = 0;
   elem(M, N, 1, 1) = 1;
}






void matrixInitZero(Matrix2x2& M)
{
   const int N = 2;

   elem(M, N, 0, 0) = 0;
   elem(M, N, 0, 1) = 0;
   elem(M, N, 1, 0) = 0;
   elem(M, N, 1, 1) = 0;
}






void matrixAdd(Matrix2x2& A, float c, const Matrix2x2& B)
{
   const int N = 2;

   elem(A, N, 0, 0) += c * elem(B, N, 0, 0);
   elem(A, N, 0, 1) += c * elem(B, N, 0, 1);
   elem(A, N, 1, 0) += c * elem(B, N, 1, 0);
   elem(A, N, 1, 1) += c * elem(B, N, 1, 1);
}






void matrixScale(Matrix2x2& A, float c)
{
   const int N = 2;

   elem(A, N, 0, 0) *= c;
   elem(A, N, 0, 1) *= c;
   elem(A, N, 1, 0) *= c;
   elem(A, N, 1, 1) *= c;
}






void matrixInverse(const Matrix2x2& A, Matrix2x2& B, float *p_det)
{
   const int N = 2;
   const float EPSILON = 1e-5;

   float det = elem(A, N, 0, 0) * elem(A, N, 1, 1) - elem(A, N, 0, 1) * elem(A, N, 1, 0);

   if ( fabs(det) <= EPSILON )
   {
      throw std::runtime_error("singular matrix");
   }

   elem(B, N, 0, 0) = +elem(A, N, 1, 1) / det;
   elem(B, N, 0, 1) = -elem(A, N, 0, 1) / det;
   elem(B, N, 1, 0) = -elem(A, N, 1, 0) / det;
   elem(B, N, 1, 1) = +elem(A, N, 0, 0) / det;

   *p_det = det;
}






void matrixProduct(const Matrix2x2& A, const Vector2& x, Vector2& b)
{
   const int N = 2;

   b.data[0] = elem(A, N, 0, 0) * x.data[0] + elem(A, N, 0, 1) * x.data[1];
   b.data[1] = elem(A, N, 1, 0) * x.data[0] + elem(A, N, 1, 1) * x.data[1];
}






void matrixOuterProduct(const Vector2& a, const Vector2& b, Matrix2x2& C)
{
   const int N = 2;

   elem(C, N, 0, 0) = a.data[0] * b.data[0];
   elem(C, N, 0, 1) = a.data[0] * b.data[1];
   elem(C, N, 1, 0) = a.data[1] * b.data[0];
   elem(C, N, 1, 1) = a.data[1] * b.data[1];
}






}
