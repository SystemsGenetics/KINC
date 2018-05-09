#include "genepair_linalg.h"



namespace Pairwise {






inline const float& elem(const Matrix2x2& M, int i, int j)
{
   return M.s[i * 2 + j];
}






inline float& elem(Matrix2x2& M, int i, int j)
{
   return M.s[i * 2 + j];
}






void vectorInitZero(Vector2& a)
{
   a.s[0] = 0;
   a.s[1] = 0;
}






void vectorAdd(Vector2& a, const Vector2& b)
{
   a.s[0] += b.s[0];
   a.s[1] += b.s[1];
}






void vectorAdd(Vector2& a, float c, const Vector2& b)
{
   a.s[0] += c * b.s[0];
   a.s[1] += c * b.s[1];
}






void vectorSubtract(Vector2& a, const Vector2& b)
{
   a.s[0] -= b.s[0];
   a.s[1] -= b.s[1];
}






void vectorScale(Vector2& a, float c)
{
   a.s[0] *= c;
   a.s[1] *= c;
}






float vectorDot(const Vector2& a, const Vector2& b)
{
   return a.s[0] * b.s[0] + a.s[1] * b.s[1];
}






float vectorDiffNorm(const Vector2& a, const Vector2& b)
{
   float dist = 0;
   dist += (a.s[0] - b.s[0]) * (a.s[0] - b.s[0]);
   dist += (a.s[1] - b.s[1]) * (a.s[1] - b.s[1]);

   return sqrt(dist);
}






void matrixInitIdentity(Matrix2x2& M)
{
   elem(M, 0, 0) = 1;
   elem(M, 0, 1) = 0;
   elem(M, 1, 0) = 0;
   elem(M, 1, 1) = 1;
}






void matrixInitZero(Matrix2x2& M)
{
   elem(M, 0, 0) = 0;
   elem(M, 0, 1) = 0;
   elem(M, 1, 0) = 0;
   elem(M, 1, 1) = 0;
}






void matrixAdd(Matrix2x2& A, float c, const Matrix2x2& B)
{
   elem(A, 0, 0) += c * elem(B, 0, 0);
   elem(A, 0, 1) += c * elem(B, 0, 1);
   elem(A, 1, 0) += c * elem(B, 1, 0);
   elem(A, 1, 1) += c * elem(B, 1, 1);
}






void matrixScale(Matrix2x2& A, float c)
{
   elem(A, 0, 0) *= c;
   elem(A, 0, 1) *= c;
   elem(A, 1, 0) *= c;
   elem(A, 1, 1) *= c;
}






void matrixInverse(const Matrix2x2& A, Matrix2x2& B, float *p_det)
{
   float det = elem(A, 0, 0) * elem(A, 1, 1) - elem(A, 0, 1) * elem(A, 1, 0);

   elem(B, 0, 0) = +elem(A, 1, 1) / det;
   elem(B, 0, 1) = -elem(A, 0, 1) / det;
   elem(B, 1, 0) = -elem(A, 1, 0) / det;
   elem(B, 1, 1) = +elem(A, 0, 0) / det;

   *p_det = det;
}






void matrixProduct(const Matrix2x2& A, const Vector2& x, Vector2& b)
{
   b.s[0] = elem(A, 0, 0) * x.s[0] + elem(A, 0, 1) * x.s[1];
   b.s[1] = elem(A, 1, 0) * x.s[0] + elem(A, 1, 1) * x.s[1];
}






void matrixOuterProduct(const Vector2& a, const Vector2& b, Matrix2x2& C)
{
   elem(C, 0, 0) = a.s[0] * b.s[0];
   elem(C, 0, 1) = a.s[0] * b.s[1];
   elem(C, 1, 0) = a.s[1] * b.s[0];
   elem(C, 1, 1) = a.s[1] * b.s[1];
}






}
