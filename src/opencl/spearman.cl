
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
 * @param x_rank
 * @param y_rank
 */
float Spearman_computeCluster(
   __global const float *x,
   __global const float *y,
   __global const char *labels,
   int sampleSize,
   char cluster,
   int minSamples,
   __global float *x_rank,
   __global float *y_rank)
{
   // extract samples in pairwise cluster
   int n = 0;

   for ( int i = 0; i < sampleSize; ++i )
   {
      if ( labels[i] == cluster )
      {
         x_rank[n] = x[i];
         y_rank[n] = y[i];
         ++n;
      }
   }

   // get power of 2 size
   int N_pow2 = nextPower2(sampleSize);

   for ( int i = n; i < N_pow2; ++i )
   {
      x_rank[i] = INFINITY;
      y_rank[i] = INFINITY;
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // compute rank of x
      bitonicSortFF(N_pow2, x_rank, y_rank);
      computeRank(x_rank, n);

      // compute rank of y
      bitonicSortFF(N_pow2, y_rank, x_rank);
      computeRank(y_rank, n);

      // compute correlation of rank arrays
      float sumx = 0;
      float sumy = 0;
      float sumx2 = 0;
      float sumy2 = 0;
      float sumxy = 0;

      for ( int i = 0; i < n; ++i )
      {
         float x_i = x_rank[i];
         float y_i = y_rank[i];

         sumx += x_i;
         sumy += y_i;
         sumx2 += x_i * x_i;
         sumy2 += y_i * y_i;
         sumxy += x_i * y_i;
      }

      result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
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
   int numPairs,
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   char clusterSize,
   __global const char *in_labels,
   int minSamples,
   __global float *work_x,
   __global float *work_y,
   __global float *out_correlations)
{
   int i = get_global_id(0);

   if ( i >= numPairs )
   {
      return;
   }

   // initialize workspace variables
   int N_pow2 = nextPower2(sampleSize);
   int2 index = in_index[i];
   __global const float *x = &expressions[index.x * sampleSize];
   __global const float *y = &expressions[index.y * sampleSize];
   __global const char *labels = &in_labels[i * sampleSize];
   __global float *x_rank = &work_x[i * N_pow2];
   __global float *y_rank = &work_y[i * N_pow2];
   __global float *correlations = &out_correlations[i * clusterSize];

   for ( char k = 0; k < clusterSize; ++k )
   {
      correlations[k] = Spearman_computeCluster(x, y, labels, sampleSize, k, minSamples, x_rank, y_rank);
   }
}
