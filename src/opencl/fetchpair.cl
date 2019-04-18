
// #include "linalg.cl"






/*!
 * Compute the initial labels for a gene pair in an expression matrix. Samples
 * with missing values and samples that fall below the expression threshold are
 * labeled as such, all other samples are labeled as cluster 0. The number of
 * clean samples is returned.
 *
 * @param globalWorkSize
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_N
 * @param out_labels
 */
__kernel void fetchPair(
   int globalWorkSize,
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   int minExpression,
   __global int *out_N,
   __global char *out_labels)
{
   int i = get_global_id(0);

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize variables
   int2 index = in_index[i];
   __global char *labels = &out_labels[i * sampleSize];
   __global int *p_N = &out_N[i];

   // index into gene expressions
   __global const float *gene1 = &expressions[index.x * sampleSize];
   __global const float *gene2 = &expressions[index.y * sampleSize];

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
         N++;
         labels[i] = 0;
      }
   }

   // save number of clean samples
   *p_N = N;
}
