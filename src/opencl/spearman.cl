
// #include "sort.cl"






int nextPower2(int n)
{
	int pow2 = 2;
	while ( pow2 < n )
	{
		pow2 *= 2;
	}

	return pow2;
}






float computeCluster(
   __global const float2 *data,
   __global const char *labels, int N,
   char cluster,
   int minSamples,
   __global float *x,
   __global float *y,
   __global int *rank)
{
   // extract samples in gene pair cluster
   int N_pow2 = nextPower2(N);
	int n = 0;

	for ( int i = 0, j = 0; i < N; ++i )
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






__kernel void computeSpearmanBlock(
   __global const float2 *in_data,
   char K,
   __global const char *in_labels, int N,
   int minSamples,
   __global float *work_x,
   __global float *work_y,
   __global int *work_rank,
   __global float *out_correlations)
{
   int i = get_global_id(0);
	int N_pow2 = nextPower2(N);

   __global const float2 *data = &in_data[i * N];
   __global const char *labels = &in_labels[i * N];
   __global float *x = &work_x[i * N_pow2];
   __global float *y = &work_y[i * N_pow2];
   __global int *rank = &work_rank[i * N_pow2];
   __global float *correlations = &out_correlations[i * K];

   for ( char k = 0; k < K; ++k )
   {
      correlations[k] = computeCluster(data, labels, N, k, minSamples, x, y, rank);
   }
}
