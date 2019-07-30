
// #include "sort.cl"






/*!
 * Remove outliers from a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only the
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples are labeled as 0, so a cluster value of 0 should be used.
 *
 * This function returns the number of clean samples remaining in the data array,
 * including samples in other clusters.
 *
 * @param x
 * @param y
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param marker
 * @param x_sorted
 * @param y_sorted
 */
int removeOutliersCluster(
   __global const float *x,
   __global const float *y,
   __global char *labels,
   int sampleSize,
   char cluster,
   char marker,
   __global float *x_sorted,
   __global float *y_sorted)
{
   // extract samples from the given cluster into separate arrays
   int n = 0;

   for ( int i = 0; i < sampleSize; i++ )
   {
      if ( labels[i] == cluster )
      {
         x_sorted[n] = x[i];
         y_sorted[n] = y[i];
         n++;
      }
   }

   // get power of 2 size
   int N_pow2 = nextPower2(sampleSize);

   for ( int i = n; i < N_pow2; ++i )
   {
      x_sorted[i] = INFINITY;
      y_sorted[i] = INFINITY;
   }

   // return if the given cluster is empty
   if ( n == 0 )
   {
      return 0;
   }

   // sort samples for each axis
   bitonicSort(x_sorted, N_pow2);
   bitonicSort(y_sorted, N_pow2);

   // compute interquartile range and thresholds for each axis
   float Q1_x = x_sorted[n * 1 / 4];
   float Q3_x = x_sorted[n * 3 / 4];
   float T_x_min = Q1_x - 1.5f * (Q3_x - Q1_x);
   float T_x_max = Q3_x + 1.5f * (Q3_x - Q1_x);

   float Q1_y = y_sorted[n * 1 / 4];
   float Q3_y = y_sorted[n * 3 / 4];
   float T_y_min = Q1_y - 1.5f * (Q3_y - Q1_y);
   float T_y_max = Q3_y + 1.5f * (Q3_y - Q1_y);

   // mark outliers
   int numSamples = 0;

   for ( int i = 0; i < sampleSize; i++ )
   {
      // mark samples in the given cluster that are outliers on either axis
      if ( labels[i] == cluster && (x[i] < T_x_min || T_x_max < x[i] || y[i] < T_y_min || T_y_max < y[i]) )
      {
         labels[i] = marker;
      }

      // count the number of remaining samples in the entire data array
      else if ( labels[i] >= 0 )
      {
         numSamples++;
      }
   }

   // return number of remaining samples
   return numSamples;
}






/*!
 * Perform outlier removal on each cluster in a parwise data array.
 *
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param in_N
 * @param in_labels
 * @param in_K
 * @param marker
 */
__kernel void removeOutliers(
   int numPairs,
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   __global int *in_N,
   __global char *in_labels,
   __global char *in_K,
   char marker,
   __global float *work_x,
   __global float *work_y)
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
   __global int *p_N = &in_N[i];
   __global char *labels = &in_labels[i * sampleSize];
   char clusterSize = in_K[i];
   __global float *x_sorted = &work_x[i * N_pow2];
   __global float *y_sorted = &work_y[i * N_pow2];

   if ( marker == -7 )
   {
      clusterSize = 1;
   }

   // do not perform post-clustering outlier removal if there is only one cluster
   if ( marker == -8 && clusterSize <= 1 )
   {
      return;
   }

   // perform outlier removal on each cluster
   int N;

   for ( char k = 0; k < clusterSize; ++k )
   {
      N = removeOutliersCluster(x, y, labels, sampleSize, k, marker, x_sorted, y_sorted);
   }

   // save number of remaining samples
   *p_N = N;
}
