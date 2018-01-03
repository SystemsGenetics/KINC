#include "genepair_linalg.h"



namespace GenePair {






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






void vector_scale(float *a, float c)
{
   a[0] *= c;
   a[1] *= c;
}






float vector_diff_norm(const float *a, const float *b)
{
   float dist = 0;
   dist += (a[0] - b[0]) * (a[0] - b[0]);
   dist += (a[1] - b[1]) * (a[1] - b[1]);

   return sqrt(dist);
}






}
