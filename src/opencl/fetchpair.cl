
// #include "linalg.cl"






/**
 * Fetch pairwise data for a pair of genes. Samples which are nan or are
 * below a threshold are excluded.
 *
 * @param expressions
 * @param size
 * @param index
 * @param minExpression
 * @param X
 * @param labels
 * @return number of rows in X
 */
int fetchPair(
   __global const float *expressions, int size,
   int2 index,
   int minExpression,
   __global Vector2 *X,
   __global char *labels)
{
   // index into gene expressions
   __global const float *gene1 = &expressions[index.x * size];
   __global const float *gene2 = &expressions[index.y * size];

   // populate X with shared expressions of gene pair
   int numSamples = 0;

   for ( int i = 0; i < size; ++i )
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
   return numSamples;
}
