
typedef float2 Vector2;






void vectorInitZero(__global Vector2 *a)
{
   a->x = 0;
   a->y = 0;
}






void vectorAdd(__global Vector2 *a, __global const Vector2 *b)
{
   a->x += b->x;
   a->y += b->y;
}






void vectorScale(__global Vector2 *a, float c)
{
   a->x *= c;
   a->y *= c;
}






float vectorDiffNorm(__global const Vector2 *a, __global const Vector2 *b)
{
   float dist = 0;
   dist += (a->x - b->x) * (a->x - b->x);
   dist += (a->y - b->y) * (a->y - b->y);

   return sqrt(dist);
}






/**
 * Implementation of the MWC64X random number generator.
 *
 * @param state
 */
uint rand(uint2 *state)
{
   enum { A = 4294883355U };
   uint x = (*state).x, c = (*state).y;  // Unpack the state
   uint res = x ^ c;                     // Calculate the result
   uint hi = mul_hi(x, A);               // Step the RNG
   x = x * A + c;
   c = hi + (x < c);
   *state = (uint2)(x, c);               // Pack the state back up

   return res;                           // Return the next result
}





/**
 * Fetch and build data matrix X for a pair of genes, skipping any expressions
 * that are missing for either gene.
 *
 * @param expressions
 * @param size
 * @param indexA
 * @param indexB
 * @param X
 * @param labels
 * @return number of rows in X
 */
int fetchData(
   __global const float *expressions, int size,
   int indexA, int indexB,
   __global Vector2 *X,
   __global int *labels)
{
   int numSamples = 0;

   indexA *= size;
   indexB *= size;

   // build data matrix from expressions and indices
   for ( int i = 0; i < size; ++i )
   {
      if ( !isnan(expressions[indexA + i]) && !isnan(expressions[indexB + i]) )
      {
         // if both expressions exist add expressions to new lists and increment
         X[numSamples] = (Vector2) (
            expressions[indexA + i],
            expressions[indexB + i]
         );
         numSamples++;

         labels[i] = 0;
      }
      else
      {
         labels[i] = -1;
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
 * @param Mu
 * @param K
 */
float computeLogLikelihood(
   __global const Vector2 *X, int N,
   __global const int *y,
   __global const Vector2 *Mu, int K)
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

         float dist = vectorDiffNorm(&X[i], &Mu[k]);

         S += dist * dist;
      }
   }

   return -S;
}






/**
 * Compute a K-means clustering model from a dataset.
 */
void computeKmeans(
   __global const Vector2 *X, int N, int K,
   __global int *y,
   float *logL,
   __global Vector2 *Mu,
   __global int *y_next)
{
   uint2 state = (get_global_id(0), get_global_id(1));

   // initialize means randomly from X
   for ( int k = 0; k < K; ++k )
   {
      int i = rand(&state) % N;
      Mu[k] = X[i];
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
            float dist = vectorDiffNorm(&X[i], &Mu[k]);

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

         vectorInitZero(&Mu[k]);

         for ( int i = 0; i < N; ++i )
         {
            if ( y[i] == k )
            {
               vectorAdd(&Mu[k], &X[i]);
               n_k++;
            }
         }

         vectorScale(&Mu[k], 1.0f / n_k);
      }
   }

   *logL = computeLogLikelihood(X, N, y, Mu, K);
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
   int minSamples, int minClusters, int maxClusters,
   __global Vector2 *work_X,
   __global int *work_y,
   __global Vector2 *work_Mu,
   __global int *work_ynext,
   __global int *result_K,
   __global int *result_labels)
{
   // initialize workspace variables
   int i = get_global_id(0);
   __global Vector2 *X = &work_X[i * size];
   __global int *y = &work_y[i * size];
   __global Vector2 *Mu = &work_Mu[i * maxClusters];
   __global int *y_next = &work_ynext[i * size];
   __global int *bestK = &result_K[i];
   __global int *bestLabels = &result_labels[i * size];

   // fetch data matrix X from expression matrix
   int numSamples = fetchData(expressions, size, pairs[i].x, pairs[i].y, X, bestLabels);

   // initialize output
   *bestK = 0;

   // make sure minimum number of related samples is reached
   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( int K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         float logL;
         computeKmeans(X, numSamples, K, y, &logL, Mu, y_next);

         // evaluate model
         float value = computeBIC(K, logL, numSamples, 2);

         // save the best model
         if ( value < bestValue )
         {
            *bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < numSamples; ++i )
            {
               if ( bestLabels[i] != -1 )
               {
                  bestLabels[i] = y[j];
                  ++j;
               }
            }
         }
      }
   }
}
