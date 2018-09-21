#include "pairwise_linalg.h"



namespace Pairwise {






/*!
 * Return the i.j element of a matrix.
 *
 * @param M
 * @param i
 * @param j
 */
inline const float& elem(const Matrix2x2& M, int i, int j)
{
   return M.s[i * 2 + j];
}






/*!
 * Return the i.j element of a matrix.
 *
 * @param M
 * @param i
 * @param j
 */
inline float& elem(Matrix2x2& M, int i, int j)
{
   return M.s[i * 2 + j];
}






/*!
 * Initialize a vector to the zero vector.
 *
 * @param a
 */
void vectorInitZero(Vector2& a)
{
   a.s[0] = 0;
   a.s[1] = 0;
}






/*!
 * Add two vectors in-place. The result is stored in a.
 *
 * @param a
 * @param b
 */
void vectorAdd(Vector2& a, const Vector2& b)
{
   a.s[0] += b.s[0];
   a.s[1] += b.s[1];
}






/*!
 * Add two vectors in-place. The vector b is scaled by a constant c, and the
 * result is stored in a.
 *
 * @param a
 * @param c
 * @param b
 */
void vectorAdd(Vector2& a, float c, const Vector2& b)
{
   a.s[0] += c * b.s[0];
   a.s[1] += c * b.s[1];
}






/*!
 * Subtract two vectors in-place. The result is stored in a.
 *
 * @param a
 * @param b
 */
void vectorSubtract(Vector2& a, const Vector2& b)
{
   a.s[0] -= b.s[0];
   a.s[1] -= b.s[1];
}






/*!
 * Scale a vector by a constant.
 *
 * @param a
 * @param c
 */
void vectorScale(Vector2& a, float c)
{
   a.s[0] *= c;
   a.s[1] *= c;
}






/*!
 * Return the dot product of two vectors.
 *
 * @param a
 * @param b
 */
float vectorDot(const Vector2& a, const Vector2& b)
{
   return a.s[0] * b.s[0] + a.s[1] * b.s[1];
}






/*!
 * Return the Euclidean distance between two vectors.
 *
 * @param a
 * @param b
 */
float vectorDiffNorm(const Vector2& a, const Vector2& b)
{
   float dist = 0;
   dist += (a.s[0] - b.s[0]) * (a.s[0] - b.s[0]);
   dist += (a.s[1] - b.s[1]) * (a.s[1] - b.s[1]);

   return sqrt(dist);
}






/*!
 * Initialize a matrix to the identity matrix.
 *
 * @param M
 */
void matrixInitIdentity(Matrix2x2& M)
{
   elem(M, 0, 0) = 1;
   elem(M, 0, 1) = 0;
   elem(M, 1, 0) = 0;
   elem(M, 1, 1) = 1;
}






/*!
 * Initialize a matrix to the zero matrix.
 *
 * @param M
 */
void matrixInitZero(Matrix2x2& M)
{
   elem(M, 0, 0) = 0;
   elem(M, 0, 1) = 0;
   elem(M, 1, 0) = 0;
   elem(M, 1, 1) = 0;
}






/*!
 * Add two matrices in place. The matrix B is scaled by a constant c, and the
 * result is stored in A.
 *
 * @param A
 * @param c
 * @param B
 */
void matrixAdd(Matrix2x2& A, float c, const Matrix2x2& B)
{
   elem(A, 0, 0) += c * elem(B, 0, 0);
   elem(A, 0, 1) += c * elem(B, 0, 1);
   elem(A, 1, 0) += c * elem(B, 1, 0);
   elem(A, 1, 1) += c * elem(B, 1, 1);
}






/*!
 * Scale a matrix by a constant.
 *
 * @param M
 * @param c
 */
void matrixScale(Matrix2x2& A, float c)
{
   elem(A, 0, 0) *= c;
   elem(A, 0, 1) *= c;
   elem(A, 1, 0) *= c;
   elem(A, 1, 1) *= c;
}






/*!
 * Compute the inverse of A and store the result in B. Additionally, the
 * determinant is returned as a pointer argument.
 *
 * @param A
 * @param B
 * @param p_det
 */
void matrixInverse(const Matrix2x2& A, Matrix2x2& B, float *p_det)
{
   float det = elem(A, 0, 0) * elem(A, 1, 1) - elem(A, 0, 1) * elem(A, 1, 0);

   elem(B, 0, 0) = +elem(A, 1, 1) / det;
   elem(B, 0, 1) = -elem(A, 0, 1) / det;
   elem(B, 1, 0) = -elem(A, 1, 0) / det;
   elem(B, 1, 1) = +elem(A, 0, 0) / det;

   *p_det = det;
}






/*!
 * Compute the matrix-vector product A * x and store the result in b.
 *
 * @param A
 * @param x
 * @param b
 */
void matrixProduct(const Matrix2x2& A, const Vector2& x, Vector2& b)
{
   b.s[0] = elem(A, 0, 0) * x.s[0] + elem(A, 0, 1) * x.s[1];
   b.s[1] = elem(A, 1, 0) * x.s[0] + elem(A, 1, 1) * x.s[1];
}






/*!
 * Compute the outer product a * b^T and store the result in C.
 *
 * @param a
 * @param b
 * @param C
 */
void matrixOuterProduct(const Vector2& a, const Vector2& b, Matrix2x2& C)
{
   elem(C, 0, 0) = a.s[0] * b.s[0];
   elem(C, 0, 1) = a.s[0] * b.s[1];
   elem(C, 1, 0) = a.s[1] * b.s[0];
   elem(C, 1, 1) = a.s[1] * b.s[1];
}






}
