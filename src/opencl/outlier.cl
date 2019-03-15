
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
 * @param data
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param marker
 * @param x_sorted
 * @param y_sorted
 */
int removeOutliersCluster(
   __global Vector2 *data,
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
         x_sorted[n] = data[i].x;
         y_sorted[n] = data[i].y;
         n++;
      }
   }

   // return if the given cluster is empty
   if ( n == 0 )
   {
      return 0;
   }

   // sort samples for each axis
   heapSort(x_sorted, n);
   heapSort(y_sorted, n);

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
      if ( labels[i] == cluster && (data[i].x < T_x_min || T_x_max < data[i].x || data[i].y < T_y_min || T_y_max < data[i].y) )
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
 * @param globalWorkSize
 * @param in_data
 * @param in_N
 * @param in_labels
 * @param sampleSize
 * @param in_K
 * @param marker
 */
__kernel void removeOutliers(
   int globalWorkSize,
   __global Vector2 *in_data,
   __global int *in_N,
   __global char *in_labels,
   int sampleSize,
   __global char *in_K,
   char marker,
   __global float *work_xy)
{
   int i = get_global_id(0);

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize workspace variables
   __global Vector2 *data = &in_data[i * sampleSize];
   __global int *p_N = &in_N[i];
   __global char *labels = &in_labels[i * sampleSize];
   char clusterSize = in_K[i];
   __global float *x_sorted = &work_xy[(2 * i + 0) * sampleSize];
   __global float *y_sorted = &work_xy[(2 * i + 1) * sampleSize];

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
      N = removeOutliersCluster(data, labels, sampleSize, k, marker, x_sorted, y_sorted);
   }

   // save number of remaining samples
   *p_N = N;
}
