
// #include "sort.cu"






/*!
 * Remove outliers from a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only the
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples are labeled as 0, so a cluster value of 0 should be used. The data
 * array should only contain samples that have a non-negative label.
 *
 * @param data
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param marker
 * @param x_sorted
 * @param y_sorted
 */
__device__
int removeOutliersCluster(
   Vector2 *data,
   char *labels,
   int sampleSize,
   char cluster,
   char marker,
   float *x_sorted,
   float *y_sorted)
{
   // extract univariate data from the given cluster
   int n = 0;

   for ( int i = 0, j = 0; i < sampleSize; i++ )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            x_sorted[n] = data[j].x;
            y_sorted[n] = data[j].y;
            n++;
         }

         j++;
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

   // remove outliers
   int numSamples = 0;

   for ( int i = 0, j = 0; i < sampleSize; i++ )
   {
      if ( labels[i] >= 0 )
      {
         // mark samples in the given cluster that are outliers on either axis
         if ( labels[i] == cluster && (data[j].x < T_x_min || T_x_max < data[j].x || data[j].y < T_y_min || T_y_max < data[j].y) )
         {
            labels[i] = marker;
         }

         // preserve all other non-outlier samples in the data array
         else
         {
            data[numSamples] = data[j];
            numSamples++;
         }

         j++;
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
__global__
void removeOutliers(
   int globalWorkSize,
   Vector2 *in_data,
   int *in_N,
   char *in_labels,
   int sampleSize,
   char *in_K,
   char marker,
   float *work_x,
   float *work_y)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize workspace variables
   Vector2 *data = &in_data[i * sampleSize];
   int *numSamples = &in_N[i];
   char *labels = &in_labels[i * sampleSize];
   char clusterSize = in_K[i];
   float *x_sorted = &work_x[i * sampleSize];
   float *y_sorted = &work_y[i * sampleSize];

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
   for ( char k = 0; k < clusterSize; ++k )
   {
      *numSamples = removeOutliersCluster(data, labels, sampleSize, k, marker, x_sorted, y_sorted);
   }
}
