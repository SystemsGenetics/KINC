
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
 * Remove outliers from a gene in a gene pair.
 *
 * @param X
 * @param N
 * @param j
 * @param labels
 * @param cluster
 * @param marker
 */
void markOutliers(
   __global const Vector2 *X, int N, int j,
   __global char *labels, char cluster,
   char marker,
   __global float *x_sorted)
{
   // compute x_sorted = X[:, j], filtered and sorted
   int n = 0;

   for ( int i = 0; i < N; i++ )
   {
      if ( labels[i] == cluster || labels[i] == marker )
      {
         x_sorted[n] = X[i].s[j];
         n++;
      }
   }

   if ( n == 0 )
   {
      return;
   }

   heapSort(x_sorted, n);

   // compute quartiles, interquartile range, upper and lower bounds
   float Q1 = x_sorted[n * 1 / 4];
   float Q3 = x_sorted[n * 3 / 4];

   float T_min = Q1 - 1.5f * (Q3 - Q1);
   float T_max = Q3 + 1.5f * (Q3 - Q1);

   // mark outliers
   for ( int i = 0; i < N; ++i )
   {
      if ( labels[i] == cluster && (X[i].s[j] < T_min || T_max < X[i].s[j]) )
      {
         labels[i] = marker;
      }
   }
}
