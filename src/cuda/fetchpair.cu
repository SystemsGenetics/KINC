
// #include "linalg.cu"






/*!
 * Extract pairwise data from an expression matrix given a pairwise index. Samples
 * with missing values and samples that fall below the expression threshold are
 * excluded. The number of extracted samples is returned.
 *
 * @param globalWorkSize
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_X
 * @param out_N
 * @param out_labels
 */
__global__
void fetchPair(
   int globalWorkSize,
   const float *expressions,
   int sampleSize,
   const int2 *in_index,
   int minExpression,
   Vector2 *out_X,
   int *out_N,
   char *out_labels)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize variables
   int2 index = in_index[i];
   Vector2 *X = &out_X[i * sampleSize];
   char *labels = &out_labels[i * sampleSize];
   int *p_N = &out_N[i];

   // index into gene expressions
   const float *gene1 = &expressions[index.x * sampleSize];
   const float *gene2 = &expressions[index.y * sampleSize];

   // populate X with shared expressions of gene pair
   int N = 0;

   for ( int i = 0; i < sampleSize; ++i )
   {
      if ( isnan(gene1[i]) || isnan(gene2[i]) )
      {
         labels[i] = -9;
      }
      else if ( gene1[i] < minExpression || gene2[i] < minExpression )
      {
         labels[i] = -6;
      }
      else
      {
         X[i] = make_float2(gene1[i], gene2[i]);
         N++;

         labels[i] = 0;
      }
   }

   // save number of clean samples
   *p_N = N;
}
