
// #include "fetchpair.cl"
// #include "linalg.cl"
// #include "outlier.cl"






typedef struct
{
   float pi;
   Vector2 mu;
   Matrix2x2 sigma;
   Matrix2x2 sigmaInv;
   float normalizer;
} Component;






void GMM_Component_initialize(
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






bool GMM_Component_prepareCovariance(__global Component *component)
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






void GMM_Component_calcLogMvNorm(
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






void GMM_kmeans(
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






void GMM_calcLogMvNorm(
   __global const Component *components, int K,
   __global const Vector2 *X, int N,
   __global float *loggamma)
{
   for ( int k = 0; k < K; ++k )
   {
      GMM_Component_calcLogMvNorm(&components[k], X, N, &loggamma[k * N]);
   }
}






void GMM_calcLogLikelihoodAndGammaNK(
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






void GMM_calcLogGammaK(
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






float GMM_calcLogGammaSum(
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






bool GMM_performMStep(
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

      bool success = GMM_Component_prepareCovariance(&components[k]);

      if ( !success )
      {
         return false;
      }
   }

   return true;
}






void GMM_calcLabels(
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






float GMM_calcEntropy(
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






/*!
 * Compute a Gaussian mixture model from a dataset.
 */
bool GMM_fit(
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

      GMM_Component_initialize(&components[k], 1.0f / K, &X[i]);
      GMM_Component_prepareCovariance(&components[k]);
   }

   // initialize means with k-means
   GMM_kmeans(components, K, X, N, MP, counts);

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
      GMM_calcLogMvNorm(components, K, X, N, loggamma);

      prevLogL = currentLogL;
      GMM_calcLogLikelihoodAndGammaNK(logpi, K, loggamma, N, &currentLogL);

      // check for convergence
      if ( fabs(currentLogL - prevLogL) < TOLERANCE )
      {
         break;
      }

      // M step
      // Let Gamma[k] = \Sum_i gamma[k, i]
      GMM_calcLogGammaK(loggamma, N, K, logGamma);

      float logGammaSum = GMM_calcLogGammaSum(logpi, K, logGamma);

      // Update parameters
      bool success = GMM_performMStep(components, K, logpi, loggamma, logGamma, logGammaSum, X, N);

      if ( !success )
      {
         return false;
      }
   }

   // save outputs
   *logL = currentLogL;
   GMM_calcLabels(loggamma, N, K, labels);
   *entropy = GMM_calcEntropy(loggamma, N, labels);

   return true;
}






typedef enum
{
   AIC,
   BIC,
   ICL
} Criterion;






/*!
 * Compute the Akaike Information Criterion of a GMM.
 */
float GMM_computeAIC(int K, int D, float logL)
{
   int p = K * (1 + D + D * D);

   return 2 * p - 2 * logL;
}






/*!
 * Compute the Bayes Information Criterion of a GMM.
 */
float GMM_computeBIC(int K, int D, float logL, int N)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL;
}






/*!
 * Compute the Integrated Completed Likelihood of a GMM.
 */
float GMM_computeICL(int K, int D, float logL, int N, float E)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL - 2 * E;
}






/*!
 * Compute a block of GMMs given a block of gene pairs.
 *
 * For each gene pair, several models are computed and the best model
 * is selected according to a criterion (BIC). The selected K and the
 * resulting sample mask for each pair is returned.
 */
__kernel void GMM_compute(
   __global const float *expressions,
   int sampleSize,
   int minSamples,
   char minClusters,
   char maxClusters,
   Criterion criterion,
   int removePreOutliers,
   int removePostOutliers,
   __global Vector2 *work_X,
   __global int *work_N,
   __global float *work_x,
   __global float *work_y,
   __global char *work_labels,
   __global Component *work_components,
   __global Vector2 *work_MP,
   __global int *work_counts,
   __global float *work_logpi,
   __global float *work_loggamma,
   __global float *work_logGamma,
   __global char *out_K,
   __global char *out_labels)
{
   int i = get_global_id(0);

   // initialize workspace variables
   __global Vector2 *X = &work_X[i * sampleSize];
   int N = work_N[i];
   __global float *x_sorted = &work_x[i * sampleSize];
   __global float *y_sorted = &work_y[i * sampleSize];
   __global char *labels = &work_labels[i * sampleSize];
   __global Component *components = &work_components[i * maxClusters];
   __global Vector2 *MP = &work_MP[i * maxClusters];
   __global int *counts = &work_counts[i * maxClusters];
   __global float *logpi = &work_logpi[i * maxClusters];
   __global float *loggamma = &work_loggamma[i * maxClusters * sampleSize];
   __global float *logGamma = &work_logGamma[i * maxClusters];
   __global char *bestK = &out_K[i];
   __global char *bestLabels = &out_labels[i * sampleSize];

   // remove pre-clustering outliers
   if ( removePreOutliers )
   {
      markOutliers(X, N, bestLabels, 0, -7, x_sorted, y_sorted);
   }

   // perform clustering only if there are enough samples
   *bestK = 0;

   if ( N >= minSamples )
   {
      float bestValue = INFINITY;

      for ( char K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         float logL;
         float entropy;

         bool success = GMM_fit(
            X, N, K,
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
         case AIC:
            value = GMM_computeAIC(K, 2, logL);
            break;
         case BIC:
            value = GMM_computeBIC(K, 2, logL, N);
            break;
         case ICL:
            value = GMM_computeICL(K, 2, logL, N, entropy);
            break;
         }

         // save the best model
         if ( value < bestValue )
         {
            *bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < N; ++i )
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

   // remove post-clustering outliers
   if ( *bestK > 1 && removePostOutliers )
   {
      for ( char k = 0; k < *bestK; ++k )
      {
         markOutliers(X, N, bestLabels, k, -8, x_sorted, y_sorted);
      }
   }
}
