
// #include "sort.cl"






/*!
 * Implementation of rand(), taken from POSIX example.
 *
 * @param state
 */
int rand(ulong *state)
{
   *state = (*state) * 1103515245 + 12345;
   return ((unsigned)((*state)/65536) % 32768);
}





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
int removeOutliers(
   __global Vector2 *data,
   __global char *labels,
   int sampleSize,
   char cluster,
   char marker,
   __global float *x_sorted,
   __global float *y_sorted)
{
   // extract univariate data from the given cluster
   int n = 0;

   for ( int i = 0, j = 0; i < sampleSize; i++ )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            x_sorted[n] = data[j].s[0];
            y_sorted[n] = data[j].s[1];
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
         if ( labels[i] == cluster && (data[j].s[0] < T_x_min || T_x_max < data[j].s[0] || data[j].s[1] < T_y_min || T_y_max < data[j].s[1]) )
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
