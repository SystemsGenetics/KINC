#include "genepair_linalg.h"



namespace GenePair {






#define ELEM(M, n, i, j) ((M)[(i) * (n) + (j)])






void vector_init_zero(float *a)
{
   a[0] = 0;
   a[1] = 0;
}






void vector_copy(float *a, const float *b)
{
   a[0] = b[0];
   a[1] = b[1];
}






void vector_add(float *a, const float *b)
{
   a[0] += b[0];
   a[1] += b[1];
}






void vector_add(float *a, float c, const float *b)
{
   a[0] += c * b[0];
   a[1] += c * b[1];
}






void vector_subtract(float *a, const float *b)
{
   a[0] -= b[0];
   a[1] -= b[1];
}






void vector_scale(float *a, float c)
{
   a[0] *= c;
   a[1] *= c;
}






float vector_dot(const float *a, const float *b)
{
   return a[0] * b[0] + a[1] * b[1];
}






float vector_diff_norm(const float *a, const float *b)
{
   float dist = 0;
   dist += (a[0] - b[0]) * (a[0] - b[0]);
   dist += (a[1] - b[1]) * (a[1] - b[1]);

   return sqrt(dist);
}






void matrix_init_identity(float *M)
{
   const int N = 2;

   ELEM(M, N, 0, 0) = 1;
   ELEM(M, N, 0, 1) = 0;
   ELEM(M, N, 1, 0) = 0;
   ELEM(M, N, 1, 1) = 1;
}






void matrix_init_zero(float *M)
{
   const int N = 2;

   ELEM(M, N, 0, 0) = 0;
   ELEM(M, N, 0, 1) = 0;
   ELEM(M, N, 1, 0) = 0;
   ELEM(M, N, 1, 1) = 0;
}






void matrix_add(float *A, float c, const float *B)
{
   const int N = 2;

   ELEM(A, N, 0, 0) += c * ELEM(B, N, 0, 0);
   ELEM(A, N, 0, 1) += c * ELEM(B, N, 0, 1);
   ELEM(A, N, 1, 0) += c * ELEM(B, N, 1, 0);
   ELEM(A, N, 1, 1) += c * ELEM(B, N, 1, 1);
}






void matrix_scale(float *A, float c)
{
   const int N = 2;

   ELEM(A, N, 0, 0) *= c;
   ELEM(A, N, 0, 1) *= c;
   ELEM(A, N, 1, 0) *= c;
   ELEM(A, N, 1, 1) *= c;
}






void matrix_inverse(const float *A, float *B, float *p_det)
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






void matrix_product(const float *A, const float *x, float *b)
{
   const int N = 2;

   b[0] = ELEM(A, N, 0, 0) * x[0] + ELEM(A, N, 0, 1) * x[1];
   b[1] = ELEM(A, N, 1, 0) * x[0] + ELEM(A, N, 1, 1) * x[1];
}






void matrix_outer_product(const float *a, const float *b, float *C)
{
   const int N = 2;

   ELEM(C, N, 0, 0) = a[0] * b[0];
   ELEM(C, N, 0, 1) = a[0] * b[1];
   ELEM(C, N, 1, 0) = a[1] * b[0];
   ELEM(C, N, 1, 1) = a[1] * b[1];
}






}
