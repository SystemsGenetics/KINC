
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






typedef struct
{
   __global Component *components;
   int K;
   float logL;
   float entropy;
   __global Vector2 *_Mu;
   __global int *_counts;
   __global float *_logpi;
   __global float *_loggamma;
   __global float *_logGamma;
   float _logGammaSum;
} GMM;






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
 * Initialize a mixture component with the given mixture weight and mean.
 *
 * @param component
 * @param pi
 * @param mu
 */
void GMM_Component_initialize(
   __global Component *component,
   float pi,
   __global const Vector2 *mu)
{
   // initialize mixture weight and mean
   component->pi = pi;
   component->mu = *mu;

   // initialize covariance to identity matrix
   matrixInitIdentity(&component->sigma);

   // initialize precision to zero matrix
   matrixInitZero(&component->sigmaInv);

   // initialize normalizer term to 0
   component->normalizer = 0;
}






/*!
 * Pre-compute the precision matrix and normalizer term for a mixture component.
 *
 * @param component
 */
bool GMM_Component_prepare(__global Component *component)
{
   const int D = 2;

   // compute precision (inverse of covariance)
   float det;
   matrixInverse(&component->sigma, &component->sigmaInv, &det);

   // return failure if matrix inverse failed
   if ( det <= 0 )
   {
      return false;
   }

   // compute normalizer term for multivariate normal distribution
   component->normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));

   // return success
   return true;
}






/*!
 * Compute the log of the probability density function of the multivariate normal
 * distribution conditioned on a single component for each point in X:
 *
 *   P(x|k) = exp(-0.5 * (x - mu)^T Sigma^-1 (x - mu)) / sqrt((2pi)^d det(Sigma))
 *
 * Therefore the log-probability is:
 *
 *   log(P(x|k)) = -0.5 * (x - mu)^T Sigma^-1 (x - mu) - 0.5 * (d * log(2pi) + log(det(Sigma)))
 *
 * @param component
 * @param X
 * @param N
 * @param logP
 */
void GMM_Component_computeLogProbNorm(
   __global const Component *component,
   __global const Vector2 *X, int N,
   __global float *logP)
{
   for (int i = 0; i < N; ++i)
   {
      // compute xm = (x - mu)
      Vector2 xm = X[i];
      vectorSubtract(&xm, &component->mu);

      // compute Sxm = Sigma^-1 xm
      Vector2 Sxm;
      matrixProduct(&component->sigmaInv, &xm, &Sxm);

      // compute xmSxm = xm^T Sigma^-1 xm
      float xmSxm = vectorDot(&xm, &Sxm);

      // compute log(P) = normalizer - 0.5 * xm^T * Sigma^-1 * xm
      logP[i] = component->normalizer - 0.5f * xmSxm;
   }
}






/*!
 * Initialize the mean of each component in the mixture model using k-means
 * clustering.
 *
 * @param gmm
 * @param X
 * @param N
 */
void GMM_initializeMeans(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   const int MAX_ITERATIONS = 20;
   const float TOLERANCE = 1e-3;
   float diff = 0;

   // initialize workspace
   __global Vector2 *Mu = gmm->_Mu;
   __global int *counts = gmm->_counts;

   for (int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t)
   {
      // compute mean and sample count for each component
      for (int k = 0; k < K; ++k)
      {
         vectorInitZero(&Mu[k]);
         counts[k] = 0;
      }

      for (int i = 0; i < N; ++i)
      {
         // determine the component mean which is nearest to x_i
         float min_dist = INFINITY;
         int min_k = 0;
         for (int k = 0; k < K; ++k)
         {
            float dist = vectorDiffNorm(&X[i], &gmm->components[k].mu);
            if (min_dist > dist)
            {
               min_dist = dist;
               min_k = k;
            }
         }

         // update mean and sample count
         vectorAdd(&Mu[min_k], &X[i]);
         ++counts[min_k];
      }

      // scale each mean by its sample count
      for (int k = 0; k < K; ++k)
      {
         vectorScale(&Mu[k], 1.0f / counts[k]);
      }

      // compute the total change of all means
      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vectorDiffNorm(&Mu[k], &gmm->components[k].mu);
      }
      diff /= K;

      // update component means
      for (int k = 0; k < K; ++k)
      {
         gmm->components[k].mu = Mu[k];
      }
   }
}






/*!
 * Perform the expectation step of the EM algorithm. In this step we compute
 * loggamma, the posterior probabilities for each component in the mixture model
 * and each sample in X, as well as logGamma, the sum of loggamma across samples,
 * and logGammaSum, another term used in the maximization step:
 *
 *   log(gamma_ki) = log(P(x_i|k)) - log(p(x))
 *
 *   log(p(x)) = a + log(sum(exp(logpi_k + logProb_ki) - a))
 *
 *   log(Gamma_k) = a + log(sum(exp(loggamma_ki) - a))
 *
 *   logGammaSum = a + log(sum(exp(logpi_k + logGamma_k) - a))
 *
 * @param gmm
 * @param X
 * @param N
 */
float GMM_computeEStep(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   // compute the log-probabilities for each component and each point in X
   __global float *logProb = gmm->_loggamma;

   for ( int k = 0; k < K; ++k )
   {
      GMM_Component_computeLogProbNorm(&gmm->components[k], X, N, &logProb[k * N]);
   }

   // compute loggamma and the log-likelihood
   float logL = 0.0;

   for (int i = 0; i < N; ++i)
   {
      // compute a = argmax(logpi_k + logProb_ki, k)
      float maxArg = -INFINITY;
      for (int k = 0; k < K; ++k)
      {
         float arg = gmm->_logpi[k] + logProb[k * N + i];
         if (maxArg < arg)
         {
            maxArg = arg;
         }
      }

      // compute logpx
      float sum = 0.0;
      for (int k = 0; k < K; ++k)
      {
         sum += exp(gmm->_logpi[k] + logProb[k * N + i] - maxArg);
      }

      float logpx = maxArg + log(sum);

      // update loggamma_ki
      for (int k = 0; k < K; ++k)
      {
         gmm->_loggamma[k * N + i] -= logpx;
      }

      // update log-likelihood
      logL += logpx;
   }

   // compute logGamma
   for (int k = 0; k < K; ++k)
   {
      __global const float *loggamma_k = &gmm->_loggamma[k * N];

      // compute a = argmax(loggamma_ki, i)
      float maxArg = -INFINITY;
      for (int i = 0; i < N; ++i)
      {
         float arg = loggamma_k[i];
         if (maxArg < arg)
         {
            maxArg = arg;
         }
      }

      // compute logGamma_k
      float sum = 0;
      for (int i = 0; i < N; ++i)
      {
         sum += exp(loggamma_k[i] - maxArg);
      }

      gmm->_logGamma[k] = maxArg + log(sum);
   }

   // compute logGammaSum
   float maxArg = -INFINITY;
   for (int k = 0; k < K; ++k)
   {
      float arg = gmm->_logpi[k] + gmm->_logGamma[k];
      if (maxArg < arg)
      {
         maxArg = arg;
      }
   }

   float sum = 0;
   for (int k = 0; k < K; ++k)
   {
      sum += exp(gmm->_logpi[k] + gmm->_logGamma[k] - maxArg);
   }

   gmm->_logGammaSum = maxArg + log(sum);

   // return log-likelihood
   return logL;
}






/*!
 * Perform the maximization step of the EM algorithm. In this step we update the
 * parameters of the the mixture model using loggamma, logGamma, and logGammaSum,
 * which are computed during the expectation step:
 *
 *   pi_k = exp(logpi_k + logGamma_k - logGammaSum)
 *
 *   mu_k = sum(exp(loggamma_ki) * x_i)) / exp(logGamma_k)
 *
 *   Sigma_k = sum(exp(loggamma_ki) * (x_i - mu_k) * (x_i - mu_k)^T) / exp(logGamma_k)
 *
 * @param gmm
 * @param X
 * @param N
 */
bool GMM_computeMStep(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   // update mixture weight
   for (int k = 0; k < K; ++k)
   {
      gmm->_logpi[k] += gmm->_logGamma[k] - gmm->_logGammaSum;

      gmm->components[k].pi = exp(gmm->_logpi[k]);
   }

   // convert loggamma to gamma
   for (int k = 0; k < K; ++k)
   {
      for (int i = 0; i < N; ++i)
      {
         gmm->_loggamma[k * N + i] = exp(gmm->_loggamma[k * N + i]);
      }
   }

   // convert logGamma to Gamma
   for (int k = 0; k < K; ++k)
   {
      gmm->_logGamma[k] = exp(gmm->_logGamma[k]);
   }

   // update remaining parameters
   for (int k = 0; k < K; ++k)
   {
      // update mean
      __global Vector2 *mu = &gmm->components[k].mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAddScaled(mu, gmm->_loggamma[k * N + i], &X[i]);
      }

      vectorScale(mu, 1.0f / gmm->_logGamma[k]);

      // update covariance matrix
      __global Matrix2x2 *sigma = &gmm->components[k].sigma;

      matrixInitZero(sigma);

      for (int i = 0; i < N; ++i)
      {
         // compute xm = (x_i - mu_k)
         Vector2 xm = X[i];
         vectorSubtract(&xm, mu);

         // compute Sigma_ki = gamma_ki * (x_i - mu_k) (x_i - mu_k)^T
         Matrix2x2 outerProduct;
         matrixOuterProduct(&xm, &xm, &outerProduct);

         matrixAddScaled(sigma, gmm->_loggamma[k * N + i], &outerProduct);
      }

      matrixScale(sigma, 1.0f / gmm->_logGamma[k]);

      // pre-compute precision matrix and normalizer term
      bool success = GMM_Component_prepare(&gmm->components[k]);

      // return failure if matrix inverse failed
      if ( !success )
      {
         return false;
      }
   }

   // return success
   return true;
}






/*!
 * Compute the cluster labels of a dataset using gamma:
 *
 *   y_i = argmax(gamma_ki, k)
 *
 * @param gamma
 * @param N
 * @param K
 * @param labels
 */
void GMM_computeLabels(
   __global const float *gamma, int N, int K,
   __global char *labels)
{
   for ( int i = 0; i < N; ++i )
   {
      // determine the value k for which gamma_ki is highest
      int max_k = -1;
      float max_gamma = -INFINITY;

      for ( int k = 0; k < K; ++k )
      {
         if ( max_gamma < gamma[k * N + i] )
         {
            max_k = k;
            max_gamma = gamma[k * N + i];
         }
      }

      // assign x_i to cluster k
      labels[i] = max_k;
   }
}






/*!
 * Compute the entropy of the mixture model for a dataset using gamma
 * and the given cluster labels:
 *
 *   E = sum(sum(z_ki * log(gamma_ki))), z_ki = (y_i == k)
 *
 * @param gamma
 * @param N
 * @param labels
 */
float GMM_computeEntropy(
   __global const float *gamma, int N,
   __global const char *labels)
{
   float E = 0;

   for ( int i = 0; i < N; ++i )
   {
      int k = labels[i];

      E += log(gamma[k * N + i]);
   }

   return E;
}






/*!
 * Fit the mixture model to a pairwise data array and compute the output cluster
 * labels for the data. The data array should only contain clean samples.
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 * @param labels
 */
bool GMM_fit(
   GMM *gmm,
   __global const Vector2 *X, int N, int K,
   __global char *labels)
{
   ulong state = 1;

   // initialize components
   gmm->K = K;

   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture weight and randomly sampled mean
      int i = rand(&state) % N;

      GMM_Component_initialize(&gmm->components[k], 1.0f / K, &X[i]);
      GMM_Component_prepare(&gmm->components[k]);
   }

   // initialize means with k-means
   GMM_initializeMeans(gmm, X, N);

   // initialize workspace
   for (int k = 0; k < K; ++k)
   {
      gmm->_logpi[k] = log(gmm->components[k].pi);
   }

   // run EM algorithm
   const int MAX_ITERATIONS = 100;
   const float TOLERANCE = 1e-8;
   float prevLogL = -INFINITY;
   float currLogL = -INFINITY;

   for ( int t = 0; t < MAX_ITERATIONS; ++t )
   {
      // perform E step
      prevLogL = currLogL;
      currLogL = GMM_computeEStep(gmm, X, N);

      // check for convergence
      if ( fabs(currLogL - prevLogL) < TOLERANCE )
      {
         break;
      }

      // perform M step
      bool success = GMM_computeMStep(gmm, X, N);

      // return failure if M-step failed (due to matrix inverse)
      if ( !success )
      {
         return false;
      }
   }

   // save outputs
   gmm->logL = currLogL;
   GMM_computeLabels(gmm->_loggamma, N, K, labels);
   gmm->entropy = GMM_computeEntropy(gmm->_loggamma, N, labels);

   return true;
}






typedef enum
{
   AIC,
   BIC,
   ICL
} Criterion;






/*!
 * Compute the Akaike Information Criterion of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 */
float GMM_computeAIC(int K, int D, float logL)
{
   int p = K * (1 + D + D * D);

   return 2 * p - 2 * logL;
}






/*!
 * Compute the Bayesian Information Criterion of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 * @param N
 */
float GMM_computeBIC(int K, int D, float logL, int N)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL;
}






/*!
 * Compute the Integrated Completed Likelihood of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 * @param N
 * @param E
 */
float GMM_computeICL(int K, int D, float logL, int N, float E)
{
   int p = K * (1 + D + D * D);

   return log((float) N) * p - 2 * logL - 2 * E;
}






/*!
 * Determine the number of clusters in a pairwise data array. Several sub-models,
 * each one having a different number of clusters, are fit to the data and the
 * sub-model with the best criterion value is selected. The data array should
 * only contain samples that have a non-negative label.
 *
 * @param globalWorkSize
 * @param sampleSize
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param out_K
 * @param out_labels
 */
__kernel void GMM_compute(
   int globalWorkSize,
   int sampleSize,
   int minSamples,
   char minClusters,
   char maxClusters,
   Criterion criterion,
   __global Vector2 *work_X,
   __global int *work_N,
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

   if ( i >= globalWorkSize )
   {
      return;
   }

   // initialize workspace variables
   __global Vector2 *data = &work_X[i * sampleSize];
   int numSamples = work_N[i];
   __global char *labels = &work_labels[i * sampleSize];
   __global Component *components = &work_components[i * maxClusters];
   __global Vector2 *Mu = &work_MP[i * maxClusters];
   __global int *counts = &work_counts[i * maxClusters];
   __global float *logpi = &work_logpi[i * maxClusters];
   __global float *loggamma = &work_loggamma[i * maxClusters * sampleSize];
   __global float *logGamma = &work_logGamma[i * maxClusters];
   __global char *bestK = &out_K[i];
   __global char *bestLabels = &out_labels[i * sampleSize];

   // initialize GMM struct
   GMM gmm = {
      .components = components,
      ._Mu = Mu,
      ._counts = counts,
      ._logpi = logpi,
      ._loggamma = loggamma,
      ._logGamma = logGamma
   };

   // perform clustering only if there are enough samples
   *bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( char K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering sub-model
         bool success = GMM_fit(&gmm, data, numSamples, K, labels);

         if ( !success )
         {
            continue;
         }

         // compute the criterion value of the sub-model
         float value = INFINITY;

         switch (criterion)
         {
         case AIC:
            value = GMM_computeAIC(K, 2, gmm.logL);
            break;
         case BIC:
            value = GMM_computeBIC(K, 2, gmm.logL, numSamples);
            break;
         case ICL:
            value = GMM_computeICL(K, 2, gmm.logL, numSamples, gmm.entropy);
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
}
