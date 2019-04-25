
// #include "sort.cu"






/*!
 * Compute the next power of 2 which occurs after a number.
 *
 * @param n
 */
__device__
int nextPower2(int n)
{
   int pow2 = 2;
   while ( pow2 < n )
   {
      pow2 *= 2;
   }

   return pow2;
}






/*!
 * Compute the Spearman correlation of a cluster in a pairwise data array.
 *
 * @param x
 * @param y
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param minSamples
 * @param x_sorted
 * @param y_sorted
 * @param rank
 */
__device__
float Spearman_computeCluster(
   const float *x,
   const float *y,
   const char *labels,
   int sampleSize,
   char cluster,
   int minSamples,
   float *x_sorted,
   float *y_sorted,
   int *rank)
{
   // extract samples in pairwise cluster
   int N_pow2 = nextPower2(sampleSize);
   int n = 0;

   for ( int i = 0; i < sampleSize; ++i )
   {
      if ( labels[i] == cluster )
      {
         x_sorted[n] = x[i];
         y_sorted[n] = y[i];
         rank[n] = n + 1;
         ++n;
      }
   }

   for ( int i = n; i < N_pow2; ++i )
   {
      x_sorted[i] = INFINITY;
      y_sorted[i] = INFINITY;
      rank[i] = 0;
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // get new power of 2 floor size
      int n_pow2 = nextPower2(n);

      // execute two sorts that are the beginning of the spearman algorithm
      bitonicSortFF(n_pow2, x_sorted, y_sorted);
      bitonicSortFI(n_pow2, y_sorted, rank);

      // go through spearman sorted rank list and calculate difference from 1,2,3,... list
      int diff = 0;

      for ( int i = 0; i < n; ++i )
      {
         int tmp = (i + 1) - rank[i];
         diff += tmp*tmp;
      }

      // compute spearman coefficient
      result = 1.0 - 6.0 * diff / (n * (n*n - 1));
   }

   return result;
}






/*!
 * Compute the correlation of each cluster in a pairwise data array. The data array
 * should only contain the clean samples that were extracted from the expression
 * matrix, while the labels should contain all samples.
 *
 * @param globalWorkSize
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param clusterSize
 * @param in_labels
 * @param minSamples
 * @param out_correlations
 */
__global__
void Spearman_compute(
   int globalWorkSize,
   const float *expressions,
   int sampleSize,
   const int2 *in_index,
   char clusterSize,
   const char *in_labels,
   int minSamples,
   float *work_xy,
   int *work_rank,
   float *out_correlations)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize workspace variables
   int N_pow2 = nextPower2(sampleSize);
   int2 index = in_index[i];
   const float *x = &expressions[index.x * sampleSize];
   const float *y = &expressions[index.y * sampleSize];
   const char *labels = &in_labels[i * sampleSize];
   float *x_sorted = &work_xy[(2 * i + 0) * N_pow2];
   float *y_sorted = &work_xy[(2 * i + 1) * N_pow2];
   int *rank = &work_rank[i * N_pow2];
   float *correlations = &out_correlations[i * clusterSize];

   for ( char k = 0; k < clusterSize; ++k )
   {
      correlations[k] = Spearman_computeCluster(x, y, labels, sampleSize, k, minSamples, x_sorted, y_sorted, rank);
   }
}
