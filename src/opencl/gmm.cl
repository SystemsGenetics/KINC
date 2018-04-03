
typedef union
{
   float s[2];
   float2 v2;
} Vector2;

typedef union
{
   float s[4];
   float4 v4;
} Matrix2x2;






#define ELEM(M, i, j) ((M)->s[(i) * 2 + (j)])






#define vectorInitZero(a) \
   (a)->s[0] = 0; \
   (a)->s[1] = 0;






#define vectorAdd(a, b) \
   (a)->s[0] += (b)->s[0]; \
   (a)->s[1] += (b)->s[1];






#define vectorAddScaled(a, c, b) \
   (a)->s[0] += (c) * (b)->s[0]; \
   (a)->s[1] += (c) * (b)->s[1];






#define vectorSubtract(a, b) \
   (a)->s[0] -= (b)->s[0]; \
   (a)->s[1] -= (b)->s[1];






#define vectorScale(a, c) \
   (a)->s[0] *= (c); \
   (a)->s[1] *= (c);






#define vectorDot(a, b) \
   ((a)->s[0] * (b)->s[0] + (a)->s[1] * (b)->s[1])






#define SQR(x) ((x)*(x))
#define vectorDiffNorm(a, b) \
   sqrt(SQR((a)->s[0] - (b)->s[0]) + SQR((a)->s[1] - (b)->s[1]))






#define matrixInitIdentity(M) \
   ELEM(M, 0, 0) = 1; \
   ELEM(M, 0, 1) = 0; \
   ELEM(M, 1, 0) = 0; \
   ELEM(M, 1, 1) = 1;






#define matrixInitZero(M) \
   ELEM(M, 0, 0) = 0; \
   ELEM(M, 0, 1) = 0; \
   ELEM(M, 1, 0) = 0; \
   ELEM(M, 1, 1) = 0;






#define matrixAddScaled(A, c, B) \
   ELEM(A, 0, 0) += (c) * ELEM(B, 0, 0); \
   ELEM(A, 0, 1) += (c) * ELEM(B, 0, 1); \
   ELEM(A, 1, 0) += (c) * ELEM(B, 1, 0); \
   ELEM(A, 1, 1) += (c) * ELEM(B, 1, 1);






#define matrixScale(A, c) \
   ELEM(A, 0, 0) *= (c); \
   ELEM(A, 0, 1) *= (c); \
   ELEM(A, 1, 0) *= (c); \
   ELEM(A, 1, 1) *= (c);






#define matrixInverse(A, B, det) \
   *det = ELEM(A, 0, 0) * ELEM(A, 1, 1) - ELEM(A, 0, 1) * ELEM(A, 1, 0); \
   ELEM(B, 0, 0) = +ELEM(A, 1, 1) / (*det); \
   ELEM(B, 0, 1) = -ELEM(A, 0, 1) / (*det); \
   ELEM(B, 1, 0) = -ELEM(A, 1, 0) / (*det); \
   ELEM(B, 1, 1) = +ELEM(A, 0, 0) / (*det);






#define matrixProduct(A, x, b) \
   (b)->s[0] = ELEM(A, 0, 0) * (x)->s[0] + ELEM(A, 0, 1) * (x)->s[1]; \
   (b)->s[1] = ELEM(A, 1, 0) * (x)->s[0] + ELEM(A, 1, 1) * (x)->s[1];






#define matrixOuterProduct(a, b, C) \
   ELEM(C, 0, 0) = (a)->s[0] * (b)->s[0]; \
   ELEM(C, 0, 1) = (a)->s[0] * (b)->s[1]; \
   ELEM(C, 1, 0) = (a)->s[1] * (b)->s[0]; \
   ELEM(C, 1, 1) = (a)->s[1] * (b)->s[1];






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
 * @param vector
 * @param minExpression
 * @param X
 * @param labels
 * @return number of rows in X
 */
int fetchData(
   __global const float *expressions, int size,
   int2 vector,
   int minExpression,
   __global Vector2 *X,
   __global char *labels)
{
   // index into gene expressions
   __global const float *gene1 = &expressions[vector.x * size];
   __global const float *gene2 = &expressions[vector.y * size];

   // populate X with shared expressions of gene pair
   int numSamples = 0;

   for ( int i = 0; i < size; ++i )
   {
      if ( isnan(gene1[i]) || isnan(gene2[i]) )
      {
         labels[i] = -9;
      }
      else if ( gene1[i] < minExpression || gene2[i] < minExpression )
      {
         labels[i] = -6;
      }
      else
      {
         X[numSamples].v2 = (float2) ( gene1[i], gene2[i] );
         numSamples++;

         labels[i] = 0;
      }
   }

   // return size of X
   return numSamples;
}






void swap(__global float *a, __global float *b)
{
   float c = *a;
   *a = *b;
   *b = c;
}






void siftDown(__global float *array, int start, int end)
{
   int root = start;

   while ( 2 * root + 1 <= end )
   {
      int child = 2 * root + 1;
      int swp = root;

      if ( array[swp] < array[child] )
      {
         swp = child;
      }

      if ( child + 1 <= end && array[swp] < array[child + 1] )
      {
         swp = child + 1;
      }

      if ( swp == root )
      {
         return;
      }
      else
      {
         swap(&array[root], &array[swp]);
         root = swp;
      }
   }
}






void heapify(__global float *array, int n)
{
   int start = ((n-1) - 1) / 2;

   while ( start >= 0 )
   {
      siftDown(array, start, n - 1);
      start -= 1;
   }
}






/**
 * Sort an array using heapsort.
 *
 * @param array
 * @param n
 */
void sort(__global float *array, int n)
{
   heapify(array, n);

   int end = n - 1;
   while ( end > 0 )
   {
      swap(&array[end], &array[0]);
      end -= 1;

      siftDown(array, 0, end);
   }
}






/**
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

   sort(x_sorted, n);

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






typedef struct
{
   float pi;
   Vector2 mu;
   Matrix2x2 sigma;
   Matrix2x2 sigmaInv;
   float normalizer;
} Component;






void Component_initialize(
   __global Component *component,
   float pi,
   __global const Vector2 *mu)
{
   // initialize pi and mu as given
   component->pi = pi;
   component->mu = *mu;

   // Use identity covariance- assume dimensions are independent
   matrixInitIdentity(&component->sigma);

   // Initialize zero artifacts
   matrixInitZero(&component->sigmaInv);

   component->normalizer = 0;
}






bool Component_prepareCovariance(__global Component *component)
{
   const int D = 2;

   // Compute inverse of Sigma once each iteration instead of
   // repeatedly for each calcLogMvNorm execution.
   float det;
   matrixInverse(&component->sigma, &component->sigmaInv, &det);

   if ( fabs(det) <= 0 )
   {
      return false;
   }

   // Compute normalizer for multivariate normal distribution
   component->normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));

   return true;
}






void Component_calcLogMvNorm(
   __global const Component *component,
   __global const Vector2 *X, int N,
   __global float *logP)
{
   // Here we are computing the probability density function of the multivariate
   // normal distribution conditioned on a single component for the set of points
   // given by X.
   //
   // P(x|k) = exp{ -0.5 * (x - mu)^T Sigma^{-} (x - mu) } / sqrt{ (2pi)^d det(Sigma) }

   for (int i = 0; i < N; ++i)
   {
      // Let xm = (x - mu)
      Vector2 xm = X[i];
      vectorSubtract(&xm, &component->mu);

      // Compute xm^T Sxm = xm^T S^-1 xm
      Vector2 Sxm;
      matrixProduct(&component->sigmaInv, &xm, &Sxm);

      float xmSxm = vectorDot(&xm, &Sxm);

      // Compute log(P) = normalizer - 0.5 * xm^T * S^-1 * xm
      logP[i] = component->normalizer - 0.5f * xmSxm;
   }
}






void kmeans(
   __global Component *components, int K,
   __global const Vector2 *X, int N,
   __global Vector2 *MP,
   __global int *counts)
{
   const int MAX_ITERATIONS = 20;
   const float TOLERANCE = 1e-3;
   float diff = 0;

   for (int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t)
   {
      // initialize old means
      for (int k = 0; k < K; ++k)
      {
         vectorInitZero(&MP[k]);
         counts[k] = 0;
      }

      // compute new means
      for (int i = 0; i < N; ++i)
      {
         float minD = INFINITY;
         int minDk = 0;
         for (int k = 0; k < K; ++k)
         {
            float dist = vectorDiffNorm(&X[i], &components[k].mu);
            if (minD > dist)
            {
               minD = dist;
               minDk = k;
            }
         }

         vectorAdd(&MP[minDk], &X[i]);
         ++counts[minDk];
      }

      for (int k = 0; k < K; ++k)
      {
         vectorScale(&MP[k], 1.0f / counts[k]);
      }

      // check for convergence
      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vectorDiffNorm(&MP[k], &components[k].mu);
      }
      diff /= K;

      // copy new means to components
      for (int k = 0; k < K; ++k)
      {
         components[k].mu = MP[k];
      }
   }
}






void calcLogMvNorm(
   __global const Component *components, int K,
   __global const Vector2 *X, int N,
   __global float *loggamma)
{
   for ( int k = 0; k < K; ++k )
   {
      Component_calcLogMvNorm(&components[k], X, N, &loggamma[k * N]);
   }
}






void calcLogLikelihoodAndGammaNK(
   __global const float *logpi, int K,
   __global float *loggamma, int N,
   float *logL)
{
   *logL = 0.0;
   for (int i = 0; i < N; ++i)
   {
      float maxArg = -INFINITY;
      for (int k = 0; k < K; ++k)
      {
         const float logProbK = logpi[k] + loggamma[k * N + i];
         if (logProbK > maxArg)
         {
            maxArg = logProbK;
         }
      }

      float sum = 0.0;
      for (int k = 0; k < K; ++k)
      {
         const float logProbK = logpi[k] + loggamma[k * N + i];
         sum += exp(logProbK - maxArg);
      }

      const float logpx = maxArg + log(sum);
      *logL += logpx;
      for (int k = 0; k < K; ++k)
      {
         loggamma[k * N + i] += -logpx;
      }
   }
}






void calcLogGammaK(
   __global const float *loggamma, int N, int K,
   __global float *logGamma)
{
   for (int k = 0; k < K; ++k)
   {
      __global const float *loggammak = &loggamma[k * N];

      float maxArg = -INFINITY;
      for (int i = 0; i < N; ++i)
      {
         const float loggammank = loggammak[i];
         if (loggammank > maxArg)
         {
            maxArg = loggammank;
         }
      }

      float sum = 0;
      for (int i = 0; i < N; ++i)
      {
         const float loggammank = loggammak[i];
         sum += exp(loggammank - maxArg);
      }

      logGamma[k] = maxArg + log(sum);
   }
}






float calcLogGammaSum(
   __global const float *logpi, int K,
   __global const float *logGamma)
{
   float maxArg = -INFINITY;
   for (int k = 0; k < K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      if (arg > maxArg)
      {
         maxArg = arg;
      }
   }

   float sum = 0;
   for (int k = 0; k < K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      sum += exp(arg - maxArg);
   }

   return maxArg + log(sum);
}






bool performMStep(
   __global Component *components, int K,
   __global float *logpi,
   __global float *loggamma,
   __global float *logGamma,
   float logGammaSum,
   __global const Vector2 *X, int N)
{
   // update pi
   for (int k = 0; k < K; ++k)
   {
      logpi[k] += logGamma[k] - logGammaSum;

      components[k].pi = exp(logpi[k]);
   }

   // convert loggamma / logGamma to gamma / Gamma to avoid duplicate exp(x) calls
   for (int k = 0; k < K; ++k)
   {
      for (int i = 0; i < N; ++i)
      {
         const int idx = k * N + i;
         loggamma[idx] = exp(loggamma[idx]);
      }
   }

   for (int k = 0; k < K; ++k)
   {
      logGamma[k] = exp(logGamma[k]);
   }

   for (int k = 0; k < K; ++k)
   {
      // Update mu
      __global Vector2 *mu = &components[k].mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAddScaled(mu, loggamma[k * N + i], &X[i]);
      }

      vectorScale(mu, 1.0f / logGamma[k]);

      // Update sigma
      __global Matrix2x2 *sigma = &components[k].sigma;

      matrixInitZero(sigma);

      for (int i = 0; i < N; ++i)
      {
         // xm = (x - mu)
         Vector2 xm = X[i];
         vectorSubtract(&xm, mu);

         // S_i = gamma_ik * (x - mu) (x - mu)^T
         Matrix2x2 outerProduct;
         matrixOuterProduct(&xm, &xm, &outerProduct);

         matrixAddScaled(sigma, loggamma[k * N + i], &outerProduct);
      }

      matrixScale(sigma, 1.0f / logGamma[k]);

      bool success = Component_prepareCovariance(&components[k]);

      if ( !success )
      {
         return false;
      }
   }

   return true;
}






void calcLabels(
   __global const float *loggamma, int N, int K,
   __global char *labels)
{
   for ( int i = 0; i < N; ++i )
   {
      int max_k = -1;
      float max_gamma = -INFINITY;

      for ( int k = 0; k < K; ++k )
      {
         if ( max_gamma < loggamma[k * N + i] )
         {
            max_k = k;
            max_gamma = loggamma[k * N + i];
         }
      }

      labels[i] = max_k;
   }
}






float calcEntropy(
   __global const float *loggamma, int N,
   __global const char *labels)
{
   float E = 0;

   for ( int i = 0; i < N; ++i )
   {
      int k = labels[i];

      E += log(loggamma[k * N + i]);
   }

   return E;
}






/**
 * Compute a Gaussian mixture model from a dataset.
 */
bool fit(
   __global const Vector2 *X, int N, int K,
   __global char *labels,
   float *logL,
   float *entropy,
   __global Component *components,
   __global Vector2 *MP,
   __global int *counts,
   __global float *logpi,
   __global float *loggamma,
   __global float *logGamma)
{
   ulong state = 1;

   // initialize components
   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture proportion and randomly sampled mean
      int i = rand(&state) % N;

      Component_initialize(&components[k], 1.0f / K, &X[i]);
      Component_prepareCovariance(&components[k]);
   }

   // initialize means with k-means
   kmeans(components, K, X, N, MP, counts);

   // initialize workspace
   for (int k = 0; k < K; ++k)
   {
      logpi[k] = log(components[k].pi);
   }

   // run EM algorithm
   const int MAX_ITERATIONS = 100;
   const float TOLERANCE = 1e-8;
   float prevLogL = -INFINITY;
   float currentLogL = -INFINITY;

   for ( int t = 0; t < MAX_ITERATIONS; ++t )
   {
      // E step
      // compute gamma, log-likelihood
      calcLogMvNorm(components, K, X, N, loggamma);

      prevLogL = currentLogL;
      calcLogLikelihoodAndGammaNK(logpi, K, loggamma, N, &currentLogL);

      // check for convergence
      if ( fabs(currentLogL - prevLogL) < TOLERANCE )
      {
         break;
      }

      // M step
      // Let Gamma[k] = \Sum_i gamma[k, i]
      calcLogGammaK(loggamma, N, K, logGamma);

      float logGammaSum = calcLogGammaSum(logpi, K, logGamma);

      // Update parameters
      bool success = performMStep(components, K, logpi, loggamma, logGamma, logGammaSum, X, N);

      if ( !success )
      {
         return false;
      }
   }

   // save outputs
   *logL = currentLogL;
   calcLabels(loggamma, N, K, labels);
   *entropy = calcEntropy(loggamma, N, labels);

   return true;
}






typedef enum
{
   BIC,
   ICL
} Criterion;






/**
 * Compute the Bayes Information Criterion of a GMM.
 */
float computeBIC(int K, float logL, int N, int D)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL;
}






/**
 * Compute the Integrated Completed Likelihood of a GMM.
 */
float computeICL(int K, float logL, int N, int D, float E)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL - 2 * E;
}






/**
 * Compute a block of GMMs given a block of gene pairs.
 *
 * For each gene pair, several models are computed and the best model
 * is selected according to a criterion (BIC). The selected K and the
 * resulting sample mask for each pair is returned.
 */
__kernel void computeGMMBlock(
   __global const float *expressions, int size,
   __global const int2 *pairs,
   int minSamples,
   int minExpression,
   char minClusters,
   char maxClusters,
   Criterion criterion,
   int removePreOutliers,
   int removePostOutliers,
   __global Vector2 *work_X,
   __global char *work_labels,
   __global Component *work_components,
   __global Vector2 *work_MP,
   __global int *work_counts,
   __global float *work_logpi,
   __global float *work_loggamma,
   __global float *work_logGamma,
   __global char *result_K,
   __global char *result_labels)
{
   int i = get_global_id(0);

   if ( pairs[i].x == 0 && pairs[i].y == 0 )
   {
      return;
   }

   // initialize workspace variables
   __global Vector2 *X = &work_X[i * size];
   __global char *labels = &work_labels[i * size];
   __global Component *components = &work_components[i * maxClusters];
   __global Vector2 *MP = &work_MP[i * maxClusters];
   __global int *counts = &work_counts[i * maxClusters];
   __global float *logpi = &work_logpi[i * maxClusters];
   __global float *loggamma = &work_loggamma[i * maxClusters * size];
   __global float *logGamma = &work_logGamma[i * maxClusters];
   __global char *bestK = &result_K[i];
   __global char *bestLabels = &result_labels[i * size];

   // fetch data matrix X from expression matrix
   int numSamples = fetchData(expressions, size, pairs[i], minExpression, X, bestLabels);

   // remove pre-clustering outliers
   __global float *work = loggamma;

   if ( removePreOutliers )
   {
      markOutliers(X, numSamples, 0, bestLabels, 0, -7, work);
      markOutliers(X, numSamples, 1, bestLabels, 0, -7, work);
   }

   // perform clustering only if there are enough samples
   *bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( char K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         float logL;
         float entropy;

         bool success = fit(
            X, numSamples, K,
            labels, &logL, &entropy,
            components,
            MP, counts,
            logpi, loggamma, logGamma
         );

         if ( !success )
         {
            continue;
         }

         // evaluate model
         float value = INFINITY;

         switch (criterion)
         {
         case BIC:
            value = computeBIC(K, logL, numSamples, 2);
            break;
         case ICL:
            value = computeICL(K, logL, numSamples, 2, entropy);
            break;
         }

         // save the best model
         if ( value < bestValue )
         {
            *bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < numSamples; ++i )
            {
               if ( bestLabels[i] >= 0 )
               {
                  bestLabels[i] = labels[j];
                  ++j;
               }
            }
         }
      }
   }

   if ( *bestK > 1 )
   {
      // remove post-clustering outliers
      if ( removePostOutliers )
      {
         for ( char k = 0; k < *bestK; ++k )
         {
            markOutliers(X, numSamples, 0, bestLabels, k, -8, work);
            markOutliers(X, numSamples, 1, bestLabels, k, -8, work);
         }
      }
   }
}
