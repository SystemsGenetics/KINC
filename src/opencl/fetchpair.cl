
// #include "linalg.cl"






/*!
 * Extract pairwise data from an expression matrix given a pairwise index. Samples
 * with missing values and samples that fall below the expression threshold are
 * excluded. The number of extracted samples is returned.
 *
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_X
 * @param out_N
 * @param out_labels
 */
__kernel void fetchPair(
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   int minExpression,
   __global Vector2 *out_X,
   __global int *out_N,
   __global char *out_labels)
{
   int i = get_global_id(0);

   // initialize variables
   int2 index = in_index[i];
   __global Vector2 *X = &out_X[i * sampleSize];
   __global char *labels = &out_labels[i * sampleSize];
   __global int *p_numSamples = &out_N[i];

   if ( index.x == 0 && index.y == 0 )
   {
      return;
   }

   // index into gene expressions
   __global const float *gene1 = &expressions[index.x * sampleSize];
   __global const float *gene2 = &expressions[index.y * sampleSize];

   // populate X with shared expressions of gene pair
   int numSamples = 0;

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
         X[numSamples].v2 = (float2) ( gene1[i], gene2[i] );
         numSamples++;

         labels[i] = 0;
      }
   }

   // return size of X
   *p_numSamples = numSamples;
}
