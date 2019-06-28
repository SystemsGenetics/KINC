
// #include "sort.cl"






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
float Spearman_computeCluster(
   __global const float *x,
   __global const float *y,
   __global const char *labels,
   int sampleSize,
   char cluster,
   int minSamples,
   __global float *x_sorted,
   __global float *y_sorted,
   __global int *rank)
{
   // extract samples in pairwise cluster
   int n = 0;

   for ( int i = 0; i < sampleSize; ++i )
   {
      if ( labels[i] == cluster )
      {
         x_sorted[n] = x[i];
         y_sorted[n] = y[i];
         rank[n] = n;
         ++n;
      }
   }

   // get power of 2 size
   int N_pow2 = nextPower2(sampleSize);

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
      // execute two sorts that are the beginning of the spearman algorithm
      bitonicSortFF(N_pow2, x_sorted, y_sorted);
      bitonicSortFI(N_pow2, y_sorted, rank);

      // go through spearman sorted rank list and calculate difference from 1,2,3,... list
      int diff = 0;

      for ( int i = 0; i < n; ++i )
      {
         int tmp = i - rank[i];
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
__kernel void Spearman_compute(
   int globalWorkSize,
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   char clusterSize,
   __global const char *in_labels,
   int minSamples,
   __global float *work_xy,
   __global int *work_rank,
   __global float *out_correlations)
{
   int i = get_global_id(0);

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize workspace variables
   int N_pow2 = nextPower2(sampleSize);
   int2 index = in_index[i];
   __global const float *x = &expressions[index.x * sampleSize];
   __global const float *y = &expressions[index.y * sampleSize];
   __global const char *labels = &in_labels[i * sampleSize];
   __global float *x_sorted = &work_xy[(2 * i + 0) * N_pow2];
   __global float *y_sorted = &work_xy[(2 * i + 1) * N_pow2];
   __global int *rank = &work_rank[i * N_pow2];
   __global float *correlations = &out_correlations[i * clusterSize];

   for ( char k = 0; k < clusterSize; ++k )
   {
      correlations[k] = Spearman_computeCluster(x, y, labels, sampleSize, k, minSamples, x_sorted, y_sorted, rank);
   }
}
