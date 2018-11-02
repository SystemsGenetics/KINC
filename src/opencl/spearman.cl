
// #include "sort.cl"






/*!
 * Compute the next power of 2 which occurs after a number.
 *
 * @param n
 */
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
 * Compute the Spearman correlation of a cluster in a pairwise data array. The
 * data array should only contain samples that have a non-negative label.
 *
 * @param data
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param minSamples
 * @param x
 * @param y
 * @param rank
 */
float Spearman_computeCluster(
   __global const float2 *data,
   __global const char *labels,
	int sampleSize,
   char cluster,
   int minSamples,
   __global float *x,
   __global float *y,
   __global int *rank)
{
   // extract samples in pairwise cluster
   int N_pow2 = nextPower2(sampleSize);
   int n = 0;

   for ( int i = 0, j = 0; i < sampleSize; ++i )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            x[n] = data[j].x;
            y[n] = data[j].y;
            rank[n] = n + 1;
            ++n;
         }

         ++j;
      }
   }

   for ( int i = n; i < N_pow2; ++i )
   {
      x[i] = INFINITY;
      y[i] = INFINITY;
      rank[i] = 0;
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // get new power of 2 floor size
      int n_pow2 = nextPower2(n);

      // execute two bitonic sorts that is beginning of spearman algorithm
      bitonicSortFF(n_pow2, x, y);
      bitonicSortFI(n_pow2, y, rank);

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
 * @param in_data
 * @param clusterSize
 * @param in_labels
 * @param sampleSize
 * @param minSamples
 * @param out_correlations
 */
__kernel void Spearman_compute(
	int globalWorkSize,
   __global const float2 *in_data,
   char clusterSize,
   __global const char *in_labels,
	int sampleSize,
   int minSamples,
   __global float *work_x,
   __global float *work_y,
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
   __global const float2 *data = &in_data[i * sampleSize];
   __global const char *labels = &in_labels[i * sampleSize];
   __global float *x = &work_x[i * N_pow2];
   __global float *y = &work_y[i * N_pow2];
   __global int *rank = &work_rank[i * N_pow2];
   __global float *correlations = &out_correlations[i * clusterSize];

   for ( char k = 0; k < clusterSize; ++k )
   {
      correlations[k] = Spearman_computeCluster(data, labels, sampleSize, k, minSamples, x, y, rank);
   }
}
