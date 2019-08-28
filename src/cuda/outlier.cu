
// #include "sort.cu"






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
__device__
int removeOutliersCluster(
   const float *x,
   const float *y,
   char *labels,
   int sampleSize,
   char cluster,
   char marker,
   float *x_sorted,
   float *y_sorted)
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
 * @param x
 * @param y
 * @param sampleSize
 * @param numSamples
 * @param labels
 * @param clusterSize
 * @param marker
 * @param x_sorted
 * @param y_sorted
 */
__device__
int removeOutliers(
   const float *x,
   const float *y,
   int sampleSize,
   int numSamples,
   char *labels,
   char clusterSize,
   char marker,
   float *x_sorted,
   float *y_sorted)
{
   // do not perform post-clustering outlier removal if there is only one cluster
   if ( marker == -8 && clusterSize <= 1 )
   {
      return numSamples;
   }

   // perform outlier removal on each cluster
   for ( char k = 0; k < clusterSize; ++k )
   {
      numSamples = removeOutliersCluster(
         x, y,
         labels,
         sampleSize,
         k,
         marker,
         x_sorted,
         y_sorted
      );
   }

   // return number of remaining samples
   return numSamples;
}
