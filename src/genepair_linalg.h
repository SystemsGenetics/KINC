#ifndef GENEPAIR_LINALG_H
#define GENEPAIR_LINALG_H
#include <ace/core/AceCore.h>

namespace GenePair
{
   void vector_init_zero(float *a);
   void vector_copy(float *a, const float *b);
   void vector_add(float *a, const float *b);
   void vector_scale(float *a, float c);
   float vector_diff_norm(const float *a, const float *b);
}

#endif
