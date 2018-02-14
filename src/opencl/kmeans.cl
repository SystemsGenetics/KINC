
typedef union
{
   float s[2];
   float2 v2;
} Vector2;






#define vectorInitZero(a) \
   (a)->s[0] = 0; \
   (a)->s[1] = 0;






#define vectorAdd(a, b) \
   (a)->s[0] += (b)->s[0]; \
   (a)->s[1] += (b)->s[1];






#define vectorScale(a, c) \
   (a)->s[0] *= (c); \
   (a)->s[1] *= (c);






#define SQR(x) ((x)*(x))
#define vectorDiffNorm(a, b) \
   sqrt(SQR((a)->s[0] - (b)->s[0]) + SQR((a)->s[1] - (b)->s[1]))






/**
 * Implementation of rand(), taken from POSIX example.
 *
 * @param state
 */
int rand(ulong *state)
{
   *state = (*state) * 1103515245 + 12345;
   return ((unsigned)((*state)/65536) % 32768);
}





/**
 * Fetch and build data matrix X for a pair of genes, skipping any expressions
 * that are missing for either gene.
 *
 * @param expressions
 * @param size
 * @param pair
 * @param minExpression
 * @param X
 * @param labels
 * @return number of rows in X
 */
int fetchData(
   __global const float *expressions, int size,
   int2 pair,
   int minExpression,
   __global Vector2 *X,
   __global char *labels)
{
   int numSamples = 0;

   pair.x *= size;
   pair.y *= size;

   // build data matrix from expressions and indices
   for ( int i = 0; i < size; ++i )
   {
      float2 v = (float2) (
         expressions[pair.x + i],
         expressions[pair.y + i]
      );

      if ( isnan(v.x) || isnan(v.y) )
      {
         labels[i] = -9;
      }
      else if ( v.x < minExpression || v.y < minExpression )
      {
         labels[i] = -6;
      }
      else
      {
         X[numSamples].v2 = v;
         numSamples++;

         labels[i] = 0;
      }
   }

   // return size of X
   return numSamples;
}






/**
 * Compute the log-likelihood of a K-means model given data X.
 *
 * @param X
 * @param N
 * @param y
 * @param means
 * @param K
 */
float computeLogLikelihood(
   __global const Vector2 *X, int N,
   __global const char *y,
   __global const Vector2 *means, int K)
{
   // compute within-class scatter
   float S = 0;

   for ( int k = 0; k < K; ++k )
   {
      for ( int i = 0; i < N; ++i )
      {
         if ( y[i] != k )
         {
            continue;
         }

         float dist = vectorDiffNorm(&X[i], &means[k]);

         S += dist * dist;
      }
   }

   return -S;
}






/**
 * Compute a K-means clustering model from a dataset.
 */
void fit(
   __global const Vector2 *X, int N, int K,
   __global char *y,
   float *logL,
   __global Vector2 *means,
   __global char *y_next)
{
   ulong state = 1;

   // initialize means randomly from X
   for ( int k = 0; k < K; ++k )
   {
      int i = rand(&state) % N;
      means[k] = X[i];
   }

   // iterate K means until convergence
   while ( true )
   {
      // compute new labels
      for ( int i = 0; i < N; ++i )
      {
         // find k that minimizes norm(x_i - mu_k)
         int min_k = -1;
         float min_dist;

         for ( int k = 0; k < K; ++k )
         {
            float dist = vectorDiffNorm(&X[i], &means[k]);

            if ( min_k == -1 || dist < min_dist )
            {
               min_k = k;
               min_dist = dist;
            }
         }

         y_next[i] = min_k;
      }

      // check for convergence
      bool converged = true;

      for ( int i = 0; i < N; ++i )
      {
         if ( y[i] != y_next[i] )
         {
            converged = false;
            break;
         }
      }

      if ( converged )
      {
         break;
      }

      // update labels
      for ( int i = 0; i < N; ++i )
      {
         y[i] = y_next[i];
      }

      // update means
      for ( int k = 0; k < K; ++k )
      {
         // compute mu_k = mean of all x_i in cluster k
         int n_k = 0;

         vectorInitZero(&means[k]);

         for ( int i = 0; i < N; ++i )
         {
            if ( y[i] == k )
            {
               vectorAdd(&means[k], &X[i]);
               n_k++;
            }
         }

         vectorScale(&means[k], 1.0f / n_k);
      }
   }

   // save outputs
   *logL = computeLogLikelihood(X, N, y, means, K);
}






/**
 * Compute the Bayes information criterion of a K-means model.
 *
 * @param K
 * @param logL
 * @param N
 * @param D
 */
float computeBIC(int K, float logL, int N, int D)
{
   int p = K * D;

   return log((float) N) * p - 2 * logL;
}






/**
 * Compute a block of K-means models given a block of gene pairs.
 *
 * For each gene pair, several models are computed and the best model
 * is selected according to a criterion (BIC). The selected K and the
 * resulting sample mask for each pair is returned.
 */
__kernel void computeKmeansBlock(
   __global const float *expressions, int size,
   __global const int2 *pairs,
   int minSamples, int minExpression, int minClusters, int maxClusters,
   __global Vector2 *work_X,
   __global char *work_y,
   __global Vector2 *work_means,
   __global char *work_ynext,
   __global int *result_K,
   __global char *result_labels)
{
   int i = get_global_id(0);

   if ( pairs[i].x == 0 && pairs[i].y == 0 )
   {
      return;
   }

   // initialize workspace variables
   __global Vector2 *X = &work_X[i * size];
   __global char *y = &work_y[i * size];
   __global Vector2 *means = &work_means[i * maxClusters];
   __global char *y_next = &work_ynext[i * size];
   __global int *bestK = &result_K[i];
   __global char *bestLabels = &result_labels[i * size];

   // fetch data matrix X from expression matrix
   int numSamples = fetchData(expressions, size, pairs[i], minExpression, X, bestLabels);

   // perform clustering only if there are enough samples
   *bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( int K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         float logL;
         fit(X, numSamples, K, y, &logL, means, y_next);

         // evaluate model
         float value = computeBIC(K, logL, numSamples, 2);

         // save the best model
         if ( value < bestValue )
         {
            *bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < numSamples; ++i )
            {
               if ( bestLabels[i] >= 0 )
               {
                  bestLabels[i] = y[j];
                  ++j;
               }
            }
         }
      }
   }
}
