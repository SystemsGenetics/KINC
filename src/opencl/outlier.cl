
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
 * Mark outliers in a vector of pairwise data. Outliers are detected independently
 * on each axis using the Tukey method, and marked with the given marker. Only those
 * samples in the given cluster are used in outlier detection. For unclustered data,
 * all samples should be labeled as 0, so a cluster value of 0 should be used.
 *
 * @param data
 * @param N
 * @param labels
 * @param cluster
 * @param marker
 * @param x_sorted
 * @param y_sorted
 */
void markOutliers(
   __global const Vector2 *data,
   int N,
   __global char *labels,
   char cluster,
   char marker,
   __global float *x_sorted,
   __global float *y_sorted)
{
   // extract univariate data from the given cluster
   int n = 0;

   for ( int i = 0; i < N; i++ )
   {
      if ( labels[i] == cluster )
      {
         x_sorted[n] = data[i].s[0];
         y_sorted[n] = data[i].s[1];
         n++;
      }
   }

   // return if the given cluster is empty
   if ( n == 0 )
   {
      return;
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
   for ( int i = 0; i < N; ++i )
   {
      if ( labels[i] == cluster )
      {
         // mark samples that are outliers on either axis
         bool outlier_x = (data[i].s[0] < T_x_min || T_x_max < data[i].s[0]);
         bool outlier_y = (data[i].s[1] < T_y_min || T_y_max < data[i].s[1]);

         if ( outlier_x || outlier_y )
         {
            labels[i] = marker;
         }
      }
   }
}
