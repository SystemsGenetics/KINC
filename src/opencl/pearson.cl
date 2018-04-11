





float computeCluster(
   __global const float2 *data,
   __global const char *labels, int N,
   char cluster,
   int minSamples,
   __global float *x,
   __global float *y)
{
   // extract samples in gene pair cluster
   int n = 0;

   for ( int i = 0, j = 0; i < N; ++i )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            x[n] = data[j].x;
            y[n] = data[j].y;
            ++n;
         }

         ++j;
      }
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // compute intermediate sums
      float sumx = 0;
      float sumy = 0;
      float sumx2 = 0;
      float sumy2 = 0;
      float sumxy = 0;

      for ( int i = 0; i < n; ++i )
      {
         sumx += x[i];
         sumy += y[i];
         sumx2 += x[i] * x[i];
         sumy2 += y[i] * y[i];
         sumxy += x[i] * y[i];
      }

      // compute Pearson correlation coefficient
      result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }

   return result;
}






__kernel void calculatePearsonBlock(
   __global const float2 *in_data,
   char K,
   __global const char *in_labels, int N,
   int minSamples,
   __global float *work_x,
   __global float *work_y,
   __global float *out_correlations)
{
   int i = get_global_id(0);

   __global const float2 *data = &in_data[i * N];
   __global const char *labels = &in_labels[i * N];
   __global float *x = &work_x[i * N];
   __global float *y = &work_y[i * N];
   __global float *correlations = &out_correlations[i * K];

   for ( char k = 0; k < K; ++k )
   {
      correlations[k] = computeCluster(data, labels, N, k, minSamples, x, y);
   }
}
