#ifndef GENEPAIR_LINALG_H
#define GENEPAIR_LINALG_H
#include <ace/core/AceCore.h>

namespace GenePair
{
   void vector_init_zero(float *a);
   void vector_copy(float *a, const float *b);
   void vector_add(float *a, const float *b);
   void vector_add(float *a, float c, const float *b);
   void vector_subtract(float *a, const float *b);
   void vector_scale(float *a, float c);
   float vector_dot(const float *a, const float *b);
   float vector_diff_norm(const float *a, const float *b);

   void matrix_init_identity(float *M);
   void matrix_init_zero(float *M);
   void matrix_add(float *A, float c, const float *B);
   void matrix_scale(float *A, float c);
   void matrix_inverse(const float *A, float *B, float *p_det);
   void matrix_product(const float *A, const float *x, float *b);
   void matrix_outer_product(const float *a, const float *b, float *C);
}

#endif
